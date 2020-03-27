library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library("colorspace")

# site_occ <- read_csv("stand occupancy.csv")
edges <- read_csv("edges.csv")
# Get a vertices list
edge_pivot <- edges %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
vertices <- edge_pivot %>% group_by(species) %>% summarise(count = n()) %>% select(-count) #nodes; vertices
rm(edge_pivot)

# G0: positive LOR only
G0 <- (tbl_graph(nodes = vertices, edges = edges)
       %>% activate(edges)
       %>% filter(p_val < 0.05)
       %>% filter(!(from == to))
       %>% filter(lor > 0)
       %>% mutate(weight = lor)
       %>% activate(nodes)
       %>% mutate(degree = degree(., V(.)))
       %>% filter(degree > 5)
       %>% mutate(community = as.factor(group_infomap())))

G0 <- (G0 %>% activate(edges)
       %>% mutate(community = ifelse(.N()$community[from] == .N()$community[to], .N()$community[from], NA)))
eg0 <- G0 %>% activate(edges) %>% as_tibble(.)
no0 <- G0 %>% activate(nodes) %>% as_tibble(.)

lo <- layout_with_fr(G0)

G0 %>% ggraph(layout = lo) + 
  geom_edge_link(aes(colour = as.factor(community)), alpha = 0.6, show.legend = F) +
  scale_edge_colour_brewer(palette = "Set1", na.value = "grey50",) +
  geom_node_point(aes(colour = (community), size = degree), alpha = 0.8) +
  scale_colour_brewer(palette = "Set1", na.value = "grey50",) +
  ggtitle('G0') 

# G1: Both positive and negative LOR, weights OR
# make a continuous palette
newcol <- colorRampPalette(brewer.pal(9, "BrBG"))
ncols <- 100
BrBG <- newcol(ncols)#apply the function to get 100 colours
rm(newcol, ncols)

mu <- mean(edges$lor)
sigma <- sd(edges$lor)

G1 <- (tbl_graph(nodes = vertices, edges = edges)
       %>% activate(edges)
       %>% filter(p_val < 0.05)
       %>% filter(!(from == to))
       %>% mutate(weight = pnorm(lor, mean = mu, sd = sigma))
       %>% activate(nodes)
       %>% mutate(degree = degree(., V(.)))
       %>% filter(degree > 5))
       
im <- cluster_infomap(G1)
G1 <-(G1 %>% activate(nodes) 
      %>% mutate(community = as.factor(membership(im))))

eg1 <- G1 %>% activate(edges) %>% as_tibble(.)
no1 <- G1 %>% activate(nodes) %>% as_tibble(.)

G1 %>% ggraph(layout = "fr") + 
  geom_edge_link(aes(colour = lor), alpha = 1, show.legend = T) +
  scale_edge_color_gradient(low = "grey80", high = "tan4") +
  geom_node_point(aes(colour = (community), size = degree), alpha = 0.8) +
  scale_colour_brewer(palette = "Set1", na.value = "grey50",) +
  ggtitle('G1') 

# Attempt to remove small modules
G2 <- G1 %>% activate(nodes) %>% filter(as.numeric(community) <= 4)
G2 %>% ggraph(layout = "fr") + 
  geom_edge_link(aes(colour = lor), alpha = 1, show.legend = T) +
  scale_edge_color_gradient(low = "yellow", high = "blue") +
  geom_node_point(aes(colour = (community), size = degree), alpha = 0.8) +
  scale_colour_brewer(palette = "Set1", na.value = "grey50",) +
  ggtitle('G1') 
