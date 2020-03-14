# Libraries
library("RMySQL")
library(tidyverse)
library(networkD3)
library(igraph)
library(plotly)
library(ggpmisc)
library(DT)
library(gridExtra)

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
    where quadrat_count = 5 and species.species_id != 4 and community is not null
    and quadrat_size = "2x2";') 
  # NOTE: this extract includes "MG5", i.e. some MG5 communities where 
  # the team have not decided
  # on a sub-group.
  
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

JointContingency <- function(d, A, B) # quadrat/assembly data, species_name, species_name
{
  d <- d %>% rename("id" = 1)
  A <- A[1]
  B <- B[1]
  As <- d %>% filter(species_name == A)
  Bs <- d %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "id")
  # Get all the assemblies/quadrats, including ones with neither A nor B
  q <- d %>% distinct(id) 
  j2 <- (left_join(q, j1, by = "id") 
         # NOTE: column length q >= column length j1
         %>% mutate(X1Y1 = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(X1Y0 = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(X0Y1 = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(X0Y0 = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  return(s)
}

OddsRatio <-  function(s) # Joint contingency (colSums)
{
  # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  o_r <- log((s[1]*s[4])/(s[2]*s[3])) # returning LOG OR
  ci_low <- o_r - 1.96*se
  ci_high <- o_r + 1.96*se
  retval <- c(o_r, ci_low, ci_high)
  names(retval) <- c("odds_ratio", "ci_low", "ci_high")
  return(retval)
}

SfChi <-  function(jc)
{
  x <- matrix(unlist(jc), ncol = 2, nrow = 2, byrow = T)
  ifelse(chisq.test(x)$p.value < 0.05, "yes", "no")
}


### END FUNCTIONS

the_data <- GetTheData()
# Restrict analysis to species with more than 20 hits
sp_count <- the_data %>% group_by(species_id) %>% summarise(n = n())
sp_list <- sp_count %>% filter(n > 20)
# the_data with reduced species count
the_data <- the_data %>% filter(species_id %in% sp_list$species_id)


assembly_data <- the_data %>% select(assembly_id, species_name) %>% distinct()
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
rm(net, n2)

# For each species pair, get the assembly odds ratio.
aor <- tibble(
  A = edges$from,
  B = edges$to,
  jc1 = 0,
  jc2 = 0,
  jc3 = 0,
  jc4 = 0,
  or = 0,
  ci_low = 0,
  ci_high = 0)

for (i in seq_along(row.names(edges)))
{
  s <- JointContingency(assembly_data, edges$from[i], edges$to[i])
  aor[i, 3:6] <- s[1:4]
  aor[i, 7:9] <- OddsRatio(s)
} # Don't remove NAs at this stage

# Make quadrat odds ratio and join wih aor
qor <- tibble(
  A = edges$from,
  B = edges$to,
  jc1 = 0,
  jc2 = 0,
  jc3 = 0,
  jc4 = 0,
  or = 0,
  ci_low = 0,
  ci_high = 0)

for (i in seq_along(row.names(edges)))
{
  s <- JointContingency(quadrat_data, edges$from[i], edges$to[i])
  qor[i, 3:6] <- s[1:4]
  qor[i, 7:9] <- OddsRatio(s)
}

pwor <- (left_join(aor, qor, by = c("A", "B"))
         %>% filter(!is.na(or.y))
         %>% filter(!is.infinite(or.y)))

# Add Chisquare test result (p < 0.05)
sfx <- tibble( sfx = pwor$A) # Arbitrary character string
for (i in seq_along(row.names(pwor)))
{
  x <- matrix(unlist(pwor[i, 3:6]), ncol = 2, nrow = 2, byrow = T)
  sf <- SfChi(x) #sfx: significance Chisquare.test: p < 0.05
  sfx$sfx[i] <- sf
  # cat(i, sf, "\n")
}
pwor$sfx <- sfx$sfx
rm(sfx)
pwor <- pwor %>% filter(sfx == "yes")

# Join share_2x2 for later use
pwor <- left_join(pwor, edges, by = c("A" = "from", "B" = "to"))

# Predict binomial p (probability of co-occurence in 5 quadrats) from bernouilli b
# (probability of co-occurence in 1 quadrat)
pred_p <- tibble(A = pwor$A,
             B = pwor$B,
             p = 0,
             b = 0,
             pp = 0)

# Joint contingency tables abcd; ja for stands, jq for quadrats
ja <- pwor %>% select(jc1.x:jc4.x)
jq <- pwor %>% select(jc1.y:jc4.y)

for (i in seq_along(row.names(pred_p)))
{
  p <- ja[i,1]/sum(ja[i,])
  b <- jq[i,1]/sum(jq[i,])
  pp <- 1 - (1 - b)^5
  
  pred_p$p[i] <- unlist(p)
  pred_p$b[i] <- unlist(b)
  pred_p$pp[i] <- unlist(pp)
}

rm(ja, jq)

plt1 <- ggplot(pred_p, aes(p, pp)) +
  geom_point() +
  geom_abline(colour = "red") +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("observed P(co-occurrence)") +
  ylab("P(co-occurence) estimated from single trial")
plot(plt1)
rm