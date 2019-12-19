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
  # Not used yet. Need to check warnings.
{
  t <- d %>% filter(species_name %in% c(A, B))
  As <- t %>% filter(species_name == A)
  browser()
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


# A <- "Achillea_millefolium"
# B <- "Trifolium_repens" # Pass in as paramters
# 
# t <- the_data %>% filter(species_name %in% c(A, B))
# As <- t %>% filter(species_name == A)
# Bs <- t %>% filter(species_name == B)
# j1 <- full_join(As, Bs, by = "quadrat_id")
# # get all the quadrats, including ones with neither A nor B
# q <- the_data %>% distinct(quadrat_id) 
# j2 <- (left_join(q, j1, by = "quadrat_id") 
#        # NOTE: column length q >= column length j1
#        %>% mutate(AandB = !is.na(species_name.x) & !is.na(species_name.y))
#        %>% mutate(AnotB = !is.na(species_name.x) & is.na(species_name.y))
#        %>% mutate(BnotA = !is.na(species_name.y) & is.na(species_name.x))
#        %>% mutate(neither = is.na(species_name.x) & is.na(species_name.y)))
# s <- colSums(j2[,4:7])
# odds <- (s[1]*s[4])/(s[2]*s[3])

edges <- (full_join(the_data, the_data, by = "quadrat_id")
  %>% rename(from = species_name.x)
  %>% rename(to = species_name.y)
  %>% filter(from != to))

# Call OddsRatio function
odds_ratios <- (head(edges, 24) 
                %>% group_by(from, to) 
                %>% summarise(odds_ratio = OddsRatio(the_data, from, to)))

