# Libraries
library("RMySQL")
library(tidyverse)
library(networkD3)
library(igraph)
library(plotly)
library(rriskDistributions)

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

# General SQL query
query <- function(q)
{
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

# Load the database
GetTheData <-  function()
{
  # GET DATA FROM DB
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  
  
  q <- sprintf('select assembly_id, assembly_name, quadrat_count, community, quadrat_id, quadrat_size, visit_date, records_id, species.species_id, 
    species.species_name from assemblies
      join quadrats on quadrats.assembly_id = assemblies_id
      join visit_dates on quadrats.vd_id = visit_dates.vds_id
      join records on records.quadrat_id = quadrats_id
      join species on species.species_id = records.species_id
      # Two assemblies have 0 quadrat count; exclude A.capillaris_stolonifera;
      # exclude some odd assemblies with no assigned community
    where quadrat_count > 0 and species.species_id != 4 and community is not null
    and quadrat_size = "2x2";') 
  # NOTE: this extract includes "MG5", i.e. some MG5 communities where 
  # the team have not decided
  # on a sub-group.
  
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

# Likelihoods
Likelihoods <- function(d, X, Y)
{
  d <- d %>% rename(id = 1)
  Xs <- d %>% filter(species_name == X)
  Ys <- d %>% filter(species_name == Y)
  # Get all the assemblies/quadrats, including ones with neither X nor Y
  trials <- d %>% distinct(id) 

  binom_dat <- left_join(trials, Xs, by = "id") %>% left_join(Ys, by = "id")
  colnames(binom_dat) <- c("id", "X", "Y")
  binom_dat <- (binom_dat %>% mutate(X1 = as.numeric(factor(X)))
                %>% mutate(Y1 = as.numeric(factor(Y)))
                %>% replace_na(list(X1=0, Y1=0))
                %>% select(id, X1, Y1)
                %>% rename(X = X1, Y = Y1))
  binom_dat <- (binom_dat %>% mutate(X1Y1 = as.numeric((X==1) & (Y==1)))
                %>% mutate(X1Y0 = as.numeric((X==1) & (Y==0)))
                %>% mutate(X0Y1 = as.numeric((X==0) & (Y==1)))
                %>% mutate(X0Y0 = as.numeric((X==0) & (Y==0))))
  
  cs <- colSums(binom_dat)
  # Use 2P(x < -|log(or)|/lse) : 2*pnorm(-abs(log(or))/se)
  lor <- log((cs[4]*cs[7])/(cs[5]*cs[6]))
  lse <- sqrt((1/cs[4]) + (1/cs[5]) + (1/cs[6] + (1/cs[7])))
  p_val <- 2*pnorm(-abs(lor/lse))
  k <- cs[4] # hits, X1Y1 ie both together in same quadrat/assembly
  n <- unlist(summarise(trials, n=n()))
  b_means <- binom_dat %>% summarise_at(c("X", "Y", "X1Y1", "X1Y0", "X0Y1", "X0Y0"), mean, na.rm = TRUE)
  # Likelihoods for data given hypothesis; hypothesis (a)true (b)false
  # Probability of exactly k hits in n trials with binomial parameter p
  # h true: observed prob X1Y1
  Ldata_h <-  dbinom(cs[4], n, b_means$X1Y1) 
  # h false: independent X, Y
  Ldata_NOTh <- dbinom(cs[2], n, b_means$X) * dbinom(cs[3], n, b_means$Y)
  retval <- c(Ldata_h, Ldata_NOTh, p_val)
  names(retval) <- c("H", "NOT_H", "p_val")
  return(retval)
}
######
the_data <- GetTheData()
# Restrict analysis to species with more than 20 hits
sp_count <- the_data %>% group_by(species_id) %>% summarise(n = n())
sp_list <- sp_count %>% filter(n > 20)
# the_data with reduced species count
the_data <- the_data %>% filter(species_id %in% sp_list$species_id)
rm(sp_count, sp_list)

assembly_data <- the_data %>% select(assembly_id, species_name) %>% distinct()
# Get species pairs
quadrat_data <- the_data %>% select(quadrat_id, species_name)
edges <- (full_join(quadrat_data, quadrat_data, by = "quadrat_id")
          %>% rename(from = species_name.x)
          %>% rename(to = species_name.y)
          %>%  select(-quadrat_id)
          %>% group_by(from, to) 
          %>% summarise(share_2x2 = n())
          %>% filter(from != to))
edges <- edges %>% filter(share_2x2 > 20)
nodes <- distinct(edges, from)
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
n2 <- simplify(net, edge.attr.comb = "first")
edges <- as_data_frame(n2, what="edges")
rm(nodes, net, n2)

prior <- 0.1
# Get likelihoods and p_values
lpv <- tibble(L_ass = rep(0.0, length(edges[,1])),
              L_NOTass = rep(0.0, length(edges[,1])),
              L_q = rep(0.0, length(edges[,1])),
              L_NOTq = rep(0.0, length(edges[,1])),
              p_val = rep(0.0, length(edges[,1])))
for (i in 1:length(edges[,1]))
#for (i in seq_along(edges))
{
  l1 <- Likelihoods(assembly_data, edges$from[i], edges$to[i])
  l2 <- Likelihoods(quadrat_data, edges$from[i], edges$to[i])
  lpv$L_ass[i] <- l1[1]
  lpv$L_NOTass[i] <- l1[2]
  lpv$L_q[i] <- l2[1]
  lpv$L_NOTq[i] <- l2[2]
  lpv$p_val[i] <- l2[3]
}
edges <- bind_cols(edges, lpv)
edges <- edges %>% filter(p_val < 0.01)

# Pr_h1 <- L1[1]*prior/(L1[1]*prior + L1[2]*(1 - prior))
# Pr_h2 <- L2[1]*prior/(L2[1]*prior + L2[2]*(1 - prior))

# Consider the two hypotheses as exclusive:

Pr_h2_ex <- L2[1]*prior_ex/((L1[1]*prior_ex) + L2[1]*prior_ex)

# ### Plots
# x <- seq(0, 1, 0.01)
# # pars_h1 <- get.beta.par(p=c(0.05, 0.5, 0.95), q=c(0.05, Pr_h1, 0.95), plot=F, show.output = F)
# # pars_h2 <- get.beta.par(p=c(0.05, 0.5, 0.95), q=c(0.05, Pr_h2, 0.95), plot=F, show.output = F)
#  pars_h1_ex <- get.beta.par(p=c(0.05, 0.5, 0.95), q=c(0.05, Pr_h1_ex, 0.95), plot=F, show.output = F)
# 
# 
# y1 <- dbeta(x, pars_h1[1], pars_h1[2])
# y2 <- dbeta(x, pars_h2[1], pars_h2[2])
# y1_ex <- dbeta(x, pars_h1_ex[1], pars_h1_ex[2])
# beta_dist <- tibble(x = x,
#                     y1 = y1,
#                     y2 = y2,
#                     y3 = y1_ex)
#  
# # p1 <- ggplot(beta_dist, aes(x=x)) +
# #   geom_line(y = y1, colour = "purple") +
# #   geom_area(aes(x=ifelse(x>0 & x<1, x, 0), y1), fill="purple", alpha=0.5) +
# #   geom_line(y=y2, colour = "springgreen") +
# #   geom_area(aes(x=ifelse(x>0 & x<1, x, 0), y2), fill="springgreen", alpha=0.5) +
# #   coord_cartesian(xlim = c(0, 1.0), 
# #                   ylim = c(0, max(max(beta_dist$y1), max(beta_dist$y2))))
# # plot(p1)
# 
# p2 <- ggplot(beta_dist, aes(x=x, y=y3)) +
#    geom_line(colour = "sienna3") +
#    geom_area(aes(x=ifelse(x>0 & x<1, x, 0), y3), fill="sienna", alpha=0.5)
# plot(p2)


# 
# # q_theta1 <- qbeta(c(0.05, 0.95), post.shape1, post.shape2) # quartiles for CrI
# # 
# # # Credibility Intervals
# # CrI1 <- q.theta[which.max(theta),1]
# # CrI2 <- q.theta[which.min(theta),2]
