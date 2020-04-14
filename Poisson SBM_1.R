library("RMySQL")
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(DT)
library(blockmodels)

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
edges <- as_tibble(as_data_frame(G, what="edges"))
nodes <- as_tibble(as_data_frame(G, what="vertices"))
rm(x, e0, G)

edges$a <- 0
edges$b <- 0
edges$c <- 0
edges$d <-  0
edges$p_val <-  0
edges$or <- 0
edges$lor <- 0
pb <- txtProgressBar(min=0, max=length(row.names(edges)))
for (i in seq_along(row.names(edges)))
{
  e1 <- as.character(edges$from[i]) # First species name
  e2 <- as.character(edges$to[i]) # Second species name
  jc <- (site_occ %>% select(e1, e2) # Two vectors, 108 long
         %>% rename(from = 1, to = 2))
  t <- table(jc$from, jc$to)
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
e0 <- edges %>% filter(p_val < 0.05) %>% filter(lor > 0)
G0 <- graph_from_data_frame(d = e0, vertices = nodes, directed = F)

# Hairball
G1 <-  as_tbl_graph(G0)
# Remove isolated nodes.
G1 <- G1 %>% activate(nodes) %>% mutate(degree = degree(G1)) %>% filter(degree > 0)
G1 %>% activate(edges) %>% ggraph(layout = layout.fruchterman.reingold(G1)) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(shape = 21, alpha = 1) +
  ggtitle('G1 SBM Poisson') +
  theme_minimal()

M <- as_adjacency_matrix(G1, type = "both", attr = "a", edges = FALSE, names = F, sparse = F)
# Blockmodel Poisson
my_model <- BM_poisson("SBM",M )
my_model$estimate()
which.max(my_model$ICL)
# Get memberships
mmZ <- my_model$memberships[[which.max(my_model$ICL)]]$Z
sbm_comms <- apply(mmZ, 1, which.max)
# Label nodes by community
G1 <-(G1 %>% activate(nodes) 
      %>% mutate(sbm_comm = as.factor(sbm_comms)))
# Identify edges WITHIN communities
G1 <- (G1 %>% activate(edges)
       %>% mutate(sbm_comm = ifelse(.N()$sbm_comm[from] == .N()$sbm_comm[to], .N()$sbm_comm[from], NA)))
# Sort nodes by SBM community prior to looking at block matrix
G1 <- G1 %>% activate(nodes) %>% arrange(sbm_comm)
eg1 <- G1 %>% activate(edges) %>% as_tibble(.)
no1 <- G1 %>% activate(nodes) %>% as_tibble(.)

# Matrix
# This plot shows edges between communities as NA
pal <- brewer.pal(25, "Dark2")
no1 <- no1 %>% mutate(axis_colour = pal[sbm_comm])

text_x <- length(row.names(no1)) - 10
text_y <-  -2

hvlines <- no1 %>% group_by(sbm_comm) %>% summarise(n = n()) %>% mutate(cs = cumsum(n) + 0.5)

ggraph(G1, 'matrix', sort.by = NULL) + 
  geom_edge_point(aes(colour = as.factor(sbm_comm)), mirror = TRUE) +
  guides(edge_colour = guide_legend(title = "community", override.aes = list(edge_size = 4))) +
  geom_vline(xintercept = hvlines$cs, alpha = 0.5, colour = "grey") +
  geom_hline(yintercept = hvlines$cs, alpha = 0.5, colour = "grey") +
  scale_edge_colour_brewer(palette = "Dark2", na.value = "grey50") +
  scale_y_reverse(breaks = seq(1, 135, by = 1), labels = no1$name, "from") +
  scale_x_continuous(breaks = seq(1, 135, by = 1), labels = no1$name, "to") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4.5, angle = 90, colour = no1$axis_colour, 
                                   face = 'bold', hjust = 1)) +
  theme(axis.text.y = element_text(size = 4.5, colour = no1$axis_colour, face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('G1 SBM Poisson, no weights, positive LOR') +
  annotate("text", label = paste("ICL = ", as.integer(max(my_model$ICL))), 
           x = text_x, y = text_y, size = 3, colour = "black")

ggsave("SBM_Poisson_1.jpg", width = 20, height = 20, units = "cm")

