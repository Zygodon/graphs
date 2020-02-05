# Predict assembly odds ratio from quadrat odds ratio treating species hits combinations
# (X1Y1, X1Y0, X0Y1, X0Y0) as Bernouilli events and so calculating a 5-quadrat assembly
# joint probability matrix from the probability of at least one event on 5 trials: 
# 1 - (1-p)^5 giving a prediction of assembly level odds ratio based only on the quadrat 
# level associations.
# If quadrat level associations were unimportant we would expect 0 correlation.
# If assembly level correlations were unimportant we would expect perfect correlation.
# The plots suggest a strong contribution to assembly correlation from quadrat
# correlation.

# Libraries
library("RMySQL")
library(tidyverse)
library(networkD3)
library(igraph)
library(plotly)

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

OddsRatio <-  function(d, A, B) # quadrat/assembly data, species_name, species_name
{
  d <- d %>% rename("id" = 1)
  A <- A[1]
  B <- B[1]
  t <- d %>% filter(species_name %in% c(A, B))
  As <- t %>% filter(species_name == A)
  Bs <- t %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "id")
  # Get all the assemblies/quadrats, including ones with neither A nor B
  q <- d %>% distinct(id) 
  j2 <- (left_join(q, j1, by = "id") 
         # NOTE: column length q >= column length j1
         %>% mutate(AandB = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(AnotB = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(BnotA = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(neither = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  gt <- sum(s)
  joint_p <- s/gt
  # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  o_r <- log((s[1]*s[4])/(s[2]*s[3])) # returning LOG OR
  ci_low <- o_r - 1.96*se
  ci_high <- o_r + 1.96*se
  # Marginal sums of the joint probability distribution
  marg_x <- (s[1] + s[2])/(s[1] + s[2] + s[3] + s[4])
  marg_y <- (s[1] + s[3])/(s[1] + s[2] + s[3] + s[4])
  retval <- c(o_r, ci_low, ci_high, marg_x, marg_y, joint_p)
  names(retval) <- c("odds_ratio", "ci_low", "ci_high", "marg_x", "marg_y", 
                     "jp1", "jp2", "jp3", "jp4")
  return(retval)
}

p1in5 <- function(p, n=5)
{
  return(1 - (1-p)^n)
}

AssemblyORGivenQuadratOR <- function(pr)
# joint probability matrix                                     
{
  a <- rep(0.0, 4)
  for (i in 1:4) a[i] <- p1in5(pr[i])
  lor <- log((a[1]*a[4])/(a[2]*a[3]))
  return(lor)
}


######
the_data <- GetTheData()
# Restrict analysis to species with more than 20 hits
sp_count <- the_data %>% group_by(species_id) %>% summarise(n = n())
sp_list <- sp_count %>% filter(n > 20)
# the_data with reduced species count
the_data <- the_data %>% filter(species_id %in% sp_list$species_id)
rm(sp_count, sp_list)

# Select followed by distinct gives us just the occurrence of each species in each assembly,
# without multiple records of its occurrence.
# This stage was missing in the previous attempt at assembly_level odds ratio.
assembly_data <- the_data %>% select(assembly_id, species_name) %>% distinct()

quadrat_data <- the_data %>% select(quadrat_id, species_name)
# In order to use the odds ratio function, we need the species pairs
# It will be better to get the species pairs from the quadrat data: assemblies
# may have shared species that don't share a quadrat. All species sharing a 
# quadrat imply a shared assembly.
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

# For each species pair, get the assembly odds ratio.
n <- length(row.names(edges))
assembly_pwor <- tibble(
  odds_ratio = 1:n,
  ci_low = 1:n,
  ci_high = 1:n,
  marg_x = 1:n,
  marg_y = 1:n,
  jp1 = 1:n,
  jp2 = 1:n,
  jp3 = 1:n,
  jp4 = 1:n)
for (i in seq_along(row.names(edges)))
{
  assembly_pwor[i,1:9] <- OddsRatio(assembly_data, edges$from[i], edges$to[i])
} 

# Now take the same set of species pairs and make the quadrats odds ratios
quadrat_pwor <- tibble(
  odds_ratio = 1:n,
  ci_low = 1:n,
  ci_high = 1:n,
  marg_x = 1:n,
  marg_y = 1:n)
for (i in seq_along(row.names(edges)))
{
  quadrat_pwor[i,1:5] <- OddsRatio(quadrat_data, edges$from[i], edges$to[i])[1:5]
}
quadrat_pwor <- bind_cols(edges, quadrat_pwor)
pwor <- bind_cols(quadrat_pwor, assembly_pwor)
colnames(pwor) <- c("from", "to", "share_2x2", "quadrat_or", "quadrat_ci_low", "quadrat_ci_high", 
                    "qmx", "qmy", "assembly_or", "assembly_ci_low", "assembly_ci_high", "amx", "amy",
                    "jp1", "jp2", "jp3", "jp4")
rm(assembly_pwor, quadrat_pwor)
pwor <- pwor %>% filter(!is.infinite(assembly_or))

# Make tibble for estimated assemblyOR from quadrat contingency
# a <- tibble(obs = pwor$assembly_or, est = pwor$assembly_or)
a <- pwor %>% select(from, to, share_2x2, assembly_or)
a$est <- 0.0
colnames(a) <- c("A", "B", "share_2x2", "obs", "est")

for (i in seq_along(row.names(a)))
{
  a$est[i] <- AssemblyORGivenQuadratOR(c(pwor$jp1[i], pwor$jp2[i], pwor$jp3[i], pwor$jp4[i]))
} 
plt1 <- ggplot(a, aes(x=obs, y=est), colour = "blue") +
  geom_abline(colour = "green") +
  #geom_errorbar(ymin = log(expected_qor$ci_low), ymax = log(expected_qor$ci_high), size = 0.1, width = 0.1, colour = "green", alpha = 0.6) +
  geom_smooth(method = "auto", size = 0.5) +
  geom_point(aes(colour = share_2x2, text = paste(A, B, sep=","), alpha = 0.5)) +
  scale_colour_gradient(low = "sienna1", high = "black") +
  labs(x = "observed assembly log(odds ratio)", y = "predicted assembly log(odds ratio)") +
  theme_grey() + coord_cartesian(xlim = c(-2, 3), ylim = c(-2,3))
plotly::ggplotly(plt1)

