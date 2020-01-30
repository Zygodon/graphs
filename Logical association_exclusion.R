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

# Return hits, X&Y (X1Y1) groupe by quadrat or assembly
# d: assembly_data or quadrat_data
HitsXY <- function(d, X, Y)
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
  retval <- c(k, n, p_val)
  names(retval) <- c("hits", "trials", "p_val")
  return(retval)
}

QuadratCount =function(d, q_id) #Return the quadrat count of the assembly
# to which the quadrat belongs. d: the_data
{
 return(unlist(d %>% filter(quadrat_id == 31) %>% distinct(quadrat_count)))
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

# Get hits (X&Y) and associated data
xy <- tibble(qhits = rep(0.0, length(edges[,1])),
             qtrials = rep(0.0, length(edges[,1])),
             qp_value = rep(0.0, length(edges[,1])),
             ahits = rep(0.0, length(edges[,1])),
             atrials = rep(0.0, length(edges[,1])),
             ap_value = rep(0.0, length(edges[,1])))

for (i in 1:length(edges[,1]))
{
  hq <- HitsXY(quadrat_data, edges$from[i], edges$to[i])
  xy$qhits[i] <- hq[1]
  xy$qtrials[i] <- hq[2]
  xy$qp_value[i] <- hq[3]
  
  hq <- HitsXY(assembly_data, edges$from[i], edges$to[i])
  xy$ahits[i] <- hq[1]
  xy$atrials[i] <- hq[2]
  xy$ap_value[i] <- hq[3]
}

edges <- bind_cols(edges, xy)
edges <- edges %>% filter(qp_value < 0.01)

edges <- (edges %>% mutate(PrQ = qhits/qtrials) 
               %>% mutate(PrA = ahits/atrials) 
               %>% mutate(BinPrA = 1-dbinom(0, 5, PrQ)))

associatives <- edges %>% select(from, to, PrQ, PrA, EPrA) %>% filter(PrA < BinPrA)
exclusives <- edges %>% select(from, to, PrQ, PrA, EPrA) %>% filter(PrA >= BinPrA)

