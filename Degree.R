library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(blockmodels)


# MAIN

# Edges stored in file. See Stand occupancy.R for code
# Key: A, species-from; B, species-to; a,b,c,d, entries in the contingency table
# p_val, significance level, Fisher's Exact Text for association;
# or, odds ratio
# lor, logarithm of the odds ratio
edges <- read_csv("edges.csv")
# Get a vertices list
vertices <- as_tibble(unique(unlist(c(edges$A, edges$B), " ")))
# G0: the basic graph
edges <- edges %>% rename(sp_from = A, sp_to = B)

# vertices1 <- edges1 %>% pivot_longer(cols = c(A, B), values_to = "species") %>% distinct(species)
G0 <- (tbl_graph(nodes = vertices, edges = edges, directed = F)
      %>% activate(nodes) %>% mutate(deg = centrality_degree()))

no0 <- G0 %>% activate(nodes) %>% as_tibble(.)

p0 <- ggplot(no0, aes(deg)) +
  geom_histogram(aes(y = ..density..), binwidth = 5, colour = "black") +
  geom_density(colour = "red", fill = "red", alpha = 0.3) +
  labs(title ="", x = "node degree")
plot(p0)

site_occ <- read_csv("stand occupancy.csv")
species_hits <- tibble(species = colnames(site_occ[2:207]), hits = colSums(site_occ[,2:207]))
no0 <- no0 %>% left_join(species_hits, by = c("value" = "species"))

p1 <- ggplot(no0, aes(hits, deg))+
  geom_point() +
  geom_smooth(method = "loess")
plot(p1)
 
(edges %>% filter(sp_from == "Luzula_campestris" | sp_to == "Luzula_campestris")
              %>% filter(lor < 1) %>% summarise(n = n()))
                        