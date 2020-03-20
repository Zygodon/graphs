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

edges$a <- 0
edges$b <- 0
edges$c <- 0
edges$d <-  0
edges$p_val <-  0
edges$or <- 0
edges$lor <- 0
pb <- txtProgressBar(min=0, max=length(row.names(edges)))
# for(i in c(0, x, 1)) {Sys.sleep(0.5); setTxtProgressBar(pb, i)}
# Sys.sleep(1)

for (i in seq_along(row.names(edges)))
{
  e1 <- as.character(edges$A[i]) # First species name
  e2 <- as.character(edges$B[i]) # Second species name
  jc <- (site_occ %>% select(e1, e2) # Two vectors, 108 long
        %>% rename(A = 1, B = 2))
    t <- table(jc$A, jc$B)
    if((0 %in% t) | sum(dim(t)) < 4){
      edges$p_val[i] <- NA
      edges$or[i] <- NA
      edges$lor[i] <- NA
    } else {
      ft <- fisher.test(t, alternative = "two.sided")
      edges$a[i] <- t[2,2]
      edges$b[i] <- t[2,1]
      edges$c[i] <- t[1,2]
      edges$d[i] <- t[1,1]
      edges$p_val[i] <- ft$p.value
      edges$or[i] <- ft$estimate
      edges$lor[i] <- log(ft$estimate)
    }
    setTxtProgressBar(pb, i)
} # end for
close(pb)
rm(i, t, e1, e2, jc, ft, pb)
edges <- edges %>% filter(!is.na(or))

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

# write_csv(site_occ, path = "stand occupancy.csv")
# write_csv(edges, path = "edges.csv")

# test <- read_csv("stand occupancy.csv")
