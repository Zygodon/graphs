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

the_data <- GetTheData()

# Make the stand (assembly) based occupancy matrix site_occ. Start with d2:
d2 <- (the_data %>% select(assembly_id, species_name)
       %>% group_by(assembly_id, species_name)
       %>% summarise(n=n())
       %>% ungroup()
       %>% pivot_wider(names_from = species_name, values_from = n))
# At this point, d2 has the number of hits for each assembly and species.
# Replace anything numeric with 1, and any NA with 0
d3 <- (d2 %>% select(-assembly_id) %>% replace(., !is.na(.), 1)
       %>% replace(., is.na(.), 0)) # Replace NAs with 0)
# Insert a column with the stand IDs. I know that it is usual with occupancy
# matrices to have the species as rows and sites as columns but as I intend
# analysis by species (R analysis as opposed to Q) I find it easier for
# the moment to think of the species in columns.
cn <- colnames(d3)
d3 <- (mutate(d3, stand = paste("S", as.character(d2$assembly_id), sep = "_"))
      %>% select(stand, cn)) #Re order columns so stand 1st
# Rename and delete d3
site_occ <- d3
rm(d2, d3)

# Pairwise list of all species combinations (expand.grid useful)
x <- colnames(site_occ)[2:length(colnames(site_occ))] # Remove the "stand" column name
e0 <- as_tibble(expand.grid(x, x))
# e0 Includes B -> A as well as A -> B
# Make a graph and simplify it to remove multiple edges
G <- graph_from_data_frame(d = e0, vertices = x, directed = F)
G <- simplify(G, edge.attr.comb = "first")
edges <- as_tibble(as_data_frame(G, what="edges")) %>% rename(A = from, B = to)
rm(x, e0, G)

# Make provision for contingencies: cell names a,b,c,d
edges$a <- 0
edges$b <- 0
edges$c <- 0
edges$d <- 0
# Get the contingency table entries a,b,c,d
for (i in seq_along(row.names(edges)))
{
  e1 <- as.character(edges$A[i]) # First species name
  e2 <- as.character(edges$B[i]) # Second species name
  jc <- (site_occ %>% select(e1, e2) # Two vectors, 108 long
                  %>% rename(A = 1, B = 2))
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
# So need to make a vertices list from edges. Can't just use the original
# species list
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
  s <- unlist(edges[i, 3:6]) # a,b,c,d
  ft <- fisher.test(matrix(s, nrow = 2, ncol = 2, byrow = T) ,alternative = "two.sided")
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
E(G0)$weight <- edges$or
plot(G0, layout = l0, vertex.label=NA)


e1 <- edges %>% filter(p_val < 0.05)
G1 <- graph_from_data_frame(d = e1, vertices = vertices, directed = F)
E(G1)$color <- ifelse(e1$lor > 0, "green", "grey50")
E(G1)$weight <- NA #e1$lor
V(G1)$size = 1
V(G1)$color <- "tomato2"
isolated <-  which(degree(G1)==0)
G1 <-  delete.vertices(G1, isolated)
# l1 <- layout.fruchterman.reingold(G1) #, weights = e1$lor)
plot(G1, vertex.label=NA)
# plot(G1, layout = l1, vertex.label=NA)

# Concentrate on the positive associations:
e1a <-  e1 %>% filter(lor > 0) # Positive associations
G1a <- graph_from_data_frame(d = e1a, vertices = vertices, directed = F)
l1a <- layout.fruchterman.reingold(G1a, weights = e1a$lor)
isolated <-  which(degree(G1a)==0)
G1a <-  delete.vertices(G1a, isolated)
E(G1a)$weight <- e1a$or
plot(G1a, vertex.label=NA)
h <- hub_score(G1a)$vector
h_score <- tibble(species = names(h), score = h)
inc_edges <- incident(G1a,  V(G1a)["Lathyrus_montanus"], mode="all")

# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(G1a))
ecol[inc_edges] <- "orange"
vcol <- rep("grey40", vcount(G1a))
vcol[V(G1a)$name=="Lathyrus_montanus"] <- "gold"
V(G1a)$size = degree(G1a)/5
plot(im, G1a, layout = l1a, vertex.color=vcol, edge.color=ecol, vertex.label=NA, axes = F, rescale = T,
ylim=c(-0.4, 0.4), xlim=c(-0.4, 0.4), asp = 1, vertex.label=NA,)

# sink("im 2020-03-18.txt")
# print(communities(im))
# sink()

G1b <- delete_edges(G1a, which(crossing(im, G1a) == T))
plot(im, G1b, vertex.label=NA,)

M1 <- delete_vertices(G1b, unlist(im[2:length(im)]))
plot(M1, vertex.label = NA)

