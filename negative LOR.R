library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library("colorspace")

site_occ <- read_csv("stand occupancy.csv")
edges <- read_csv("edges.csv")
# Get a vertices list
edge_pivot <- edges %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
vertices <- edge_pivot %>% group_by(species) %>% summarise(count = n()) %>% select(-count) #nodes; vertices
rm(edge_pivot)

G0 <- (tbl_graph(nodes = vertices, edges = edges)
       %>% activate(edges)
       %>% filter(p_val < 0.05)
       %>% filter(!(from == to))
       %>% filter(lor < 0)
       %>% mutate(weight = lor)
       %>% activate(nodes)
       %>% mutate(degree = degree(., V(.)))
       %>% filter(degree > 5))
       # %>% activate(edges)
       # %>% filter(p_val < 0.05)
       # %>% filter(!(from == to))
       # %>% filter(lor < 0)
       # %>% mutate(weight = lor))

lo <- layout_with_fr(G0)

G0 %>% ggraph(layout = lo) + 
#G0 %>% ggraph(layout = 'fr') + 
  geom_edge_link(colour = "red", alpha = 0.6) + 
  geom_node_point(aes(size = degree), alpha = 0.8) + 
  # geom_node_text(aes(label=species), repel=TRUE, size =1.5,colour="black") +
  ggtitle('G0') 

