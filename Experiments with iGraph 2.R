library("RMySQL")
library(tidyverse)
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

SetOne <- function(x, na.rm = T) (1)

the_data <- GetTheData()

# Make the stand (assembly) based occupancy matrix d2
d2 <- (the_data %>% select(assembly_id, species_name)
       %>% group_by(assembly_id, species_name)
       %>% summarise(n=n())
       %>% ungroup()
       %>% pivot_wider(names_from = species_name, values_from = n))
# At this point, d2 has the number of hits for each assembly and species.
# Replace anything numeric with 1, and any NA with 0
d3 <- (d2 %>% select(-assembly_id) %>% replace(., !is.na(.), 1)
       %>% replace(., is.na(.), 0)) # Replace NAs with 0)
cn <- colnames(d3)
d3 <- (mutate(d3, stand = paste("S", as.character(d2$assembly_id), sep = "_"))
      %>% select(stand, cn)) #Re order columns so stand 1st
# Rename and delete d3
site_occ <- d3
rm(d2, d3)

x <- colnames(site_occ)[2:length(colnames(site_occ))] # Remove the "stand" column name
edges <- as_tibble(expand.grid(x, x))
rm(x)
edges <- (edges %>% rename(A = 1) 
          %>% rename(B = 2)
          %>% filter(A != B))

edges$a <- 0
edges$b <- 0
edges$c <- 0
edges$d <- 0
for (i in seq_along(row.names(edges)))
{
  e1 <- as.character(edges$A[i])
  e2 <- as.character(edges$B[i])
  jc <- site_occ %>% select(e1, e2) %>% rename(A = 1, B = 2)
  jc <- jc %>% mutate(a = ifelse(A == 1 & B == 1, 1, 0))
  jc <- jc %>% mutate(b = ifelse(A == 1 & B == 0, 1, 0))
  jc <- jc %>% mutate(c = ifelse(A == 0 & B == 1, 1, 0))
  jc <- jc %>% mutate(d = ifelse(A == 0 & B == 0, 1, 0))
  s <- colSums(jc[3:6])
  edges$a[i] <- s[1]
  edges$b[i] <- s[2]
  edges$c[i] <- s[3]
  edges$d[i] <- s[4]
}
rm(jc)

# The following operation reduces the number of edges, and potentially 
# eliminates some species.
# Avoid div 0 in (a*c)/(b*d)
# Avoid 0 numerator in (a*c)/(b*d)
edges <- (edges %>% filter(a > 0)
          %>% filter (b > 0) 
          %>% filter (c > 0)
          %>% filter (d > 0))      
# So need to make a vertices list from edges.
edge_pivot <- edges %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
vertices <- edge_pivot %>% group_by(species) %>% summarise(count = n()) %>% select(-count) #nodes; vertices
rm(edge_pivot)

# Use fisher.test to make p.val, odds ratio, log odds ratio and log odds ratio^2
edges$p_val <- 0
edges$or <- 0
edges$lor <- 0
edges$lor2 <- 0
for (i in seq_along(row.names(edges)))
{
  s <- unlist(edges[i, 3:6])
  ft <- fisher.test(matrix(s, nrow = 2, ncol = 2, byrow = T) ,alternative = "greater")
  edges[i, 7] <- ft$p.value
  edges[i, 8] <- ft$estimate            # or
  edges[i, 9] <- log(ft$estimate)       # lor
  edges[i, 10] <- (log(ft$estimate))^2  # lor^2
}
rm(ft)

# Check on the lor histogram
plt1 <- ggplot(edges, aes(lor))  +
  geom_histogram(aes(y = ..density..), binwidth = 0.25, colour = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(edges$lor), sd = sd(edges$lor)), colour = "green") +
  geom_vline(xintercept = 0, colour = "red") +
  xlim(-3, 6) +
  ylim(0.0, 0.4) +
  labs(title ="", x = "log(odds ratio)")
plot(plt1)
t.test(edges$lor, mu = 0) # p-value < 2.2e-16

G0 <- graph_from_data_frame(d = edges, vertices = vertices, directed = F)
l0 <- layout.fruchterman.reingold(G0, weights = edges$or)
# E(G0)$weight <- edges$or
plot(G0, layout = l0, vertex.label=NA)

e1 <- edges %>% filter(p_val < 0.05)
G1 <- graph_from_data_frame(d = e1, vertices = vertices, directed = F)
isolated <-  which(degree(G1)==0)
G1 <-  delete.vertices(G1, isolated)
l1 <- layout.fruchterman.reingold(G1, weights = e1$or)
plot(G1, layout = l1, vertex.label=NA)

