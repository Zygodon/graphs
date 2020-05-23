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

h <- seq(0, 100, 1)
d <- (145 * (1 - exp(-h/8.686))) - (0.959 * h)
m <- tibble(h = h, d = d)

p1 <- ggplot(no0, aes(hits, deg))+
  geom_point() +
  geom_smooth(method = "loess") +
  geom_abline(slope = 8.686, intercept = 16.3, colour = "red") +
  geom_abline(slope = -0.959, intercept = 186.7, colour = "red") +
  geom_smooth(data = m, aes(h, d), colour = "green")
plot(p1)
 
# Parameters need optimising but approach looks OK
# i. e. could model the fit by an exponential representing the
# inevitable increase in degree with commonness plus a negative linear 
# relationship between degree and the area occupied by very common plants.


# (edges %>% filter(sp_from == "Luzula_campestris" | sp_to == "Luzula_campestris")
#               %>% filter(lor < 1) %>% summarise(n = n()))
                        