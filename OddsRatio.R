# Libraries
library("RMySQL")
library(tidyverse)
# library(corrplot)
library(networkD3)
library(igraph)

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
  
  
  q <- sprintf('select assembly_id, assembly_name, quadrat_count, community, quadrat_id, visit_date, records_id, species.species_id, 
    species.species_name from assemblies
      join quadrats on quadrats.assembly_id = assemblies_id
      join visit_dates on quadrats.vd_id = visit_dates.vds_id
      join records on records.quadrat_id = quadrats_id
      join species on species.species_id = records.species_id
      # Two assemblies have 0 quadrat count; exclude A.capillaris_stolonifera;
      # exclude some odd assemblies with no assigned community
    where quadrat_count > 0 and species.species_id != 4 and community is not null;') 
  # NOTE: this extract includes "MG5", i.e. some MG5 communities where 
  # the team have not decided
  # on a sub-group.
  
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

OddsRatio <-  function(d, A, B) # the_data, species_name, species_name
{
  A <- A[1]
  B <- B[1]
  t <- d %>% filter(species_name %in% c(A, B))
  As <- t %>% filter(species_name == A)
  Bs <- t %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "quadrat_id")
  # get all the quadrats, including ones with neither A nor B
  q <- d %>% distinct(quadrat_id) 
  j2 <- (left_join(q, j1, by = "quadrat_id") 
         # NOTE: column length q >= column length j1
         %>% mutate(AandB = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(AnotB = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(BnotA = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(neither = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  o_r <- (s[1]*s[4])/(s[2]*s[3])
  ci_low <- exp(log(o_r)-1.96*se)
  ci_high <- exp(log(o_r)+1.96*se)
  retval <- c(o_r, ci_low, ci_high)
  names(retval) <- c("odds_ratio", "ci_low", "ci_high")
  return(retval)
}

# End of functions

d <- GetTheData()
# the_data <- d %>% filter(community == "MG5c")
the_data <- (d %>% filter(grepl("MG", community))
             %>% select(quadrat_id, species_name))
rm(d)

# Odds ratio function
OddsRatio <-  function(d, A, B) # the_data, species_name, species_name
{
  A <- A[1] # Avoiding Warnings
  B <- B[1]
  t <- d %>% filter(species_name %in% c(A, B))
  As <- t %>% filter(species_name == A)
  Bs <- t %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "quadrat_id")
  # get all the quadrats, including ones with neither A nor B
  q <- d %>% distinct(quadrat_id) 
  j2 <- (left_join(q, j1, by = "quadrat_id") 
         # NOTE: column length q >= column length j1
         %>% mutate(AandB = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(AnotB = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(BnotA = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(neither = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  return((s[1]*s[4])/(s[2]*s[3]))
}

edges <- (full_join(the_data, the_data, by = "quadrat_id")
           %>% rename(from = species_name.x)
           %>% rename(to = species_name.y)
           %>%  select(-quadrat_id)
           %>% group_by(from, to) 
           %>% summarise(wt = n())
           %>% filter(from != to))
nodes <- distinct(edges, from)
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
n2 <- simplify(net, edge.attr.comb = "ignore")
edges2 <- as_data_frame(n2, what="edges")

n <- length(row.names(edges2))

pairwise_o_r <- tibble(
  odds_ratio = 1:n,
  ci_low = 1:n,
  ci_high = 1:n
)

for (i in seq_along(row.names(edges2)))
{
  pairwise_o_r[i,1:3] <- OddsRatio(the_data, edges2$from[i], edges2$to[i])
}
edges3 <-( bind_cols(edges2, pairwise_o_r)
           %>% filter(!is.infinite(odds_ratio))
           %>% filter(!is.infinite(ci_low))
           %>% filter(!is.infinite(ci_high))
           %>% filter(!is.na(ci_low))
           %>% filter(!is.na(ci_high))
           %>% filter(ci_low > 1))
net2 <- graph_from_data_frame(d = edges3, vertices = nodes, directed = F)
plot(net2)
l <- layout.fruchterman.reingold(net2)
plot(net2, layout=l, vertex.label=NA)
