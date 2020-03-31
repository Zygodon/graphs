library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)

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
cle <- cluster_leading_eigen(G0, weights = NA)
G0 <-(G0 %>% activate(nodes) 
      %>% mutate(community = as.factor(membership(cle))))

lcle <- layout_with_fr(G0)
G0 %>% activate(edges) %>% ggraph(layout = lcle) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(aes(fill = community, size = degree), show.legend = T, shape = 21, alpha = 1) +
  scale_fill_brewer(palette = "Dark2", na.value = "grey50",) +
  ggtitle('G0 Cluster Leading Eigen, no weights, positive LOR') +
  theme_minimal() + th_no_axes( ) +
  guides(fill = guide_legend(override.aes = list(size=5))) # +

cle_clusters <- as_tibble(G0 %>% activate(nodes))




# mu <- mean(edges$lor)
# sigma <- sd(edges$lor)
# 
# G1 <- (tbl_graph(nodes = vertices, edges = edges)
#        %>% activate(edges)
#        %>% filter(p_val < 0.05)
#        %>% filter(!(from == to))
#        %>% mutate(weight = pnorm(lor, mean = mu, sd = sigma))
#        %>% activate(nodes)
#        %>% mutate(degree = degree(., V(.)))
#        %>% filter(degree > 5))

### CLUSTER INFOMAP       
# im <- cluster_infomap(G1)
# G1 <-(G1 %>% activate(nodes) 
#       %>% mutate(community = as.factor(membership(im))))
# 
# # Remove small modules
# G2 <- G1 %>% activate(nodes) %>% filter(as.numeric(community) <= 3)
# eg2 <- G2 %>% activate(edges) %>% as_tibble(.)
# no2 <- G2 %>% activate(nodes) %>% as_tibble(.)
# 
# lg2 <- layout_with_fr(G2)
# 
# plt <- ggraph(G2, layout = lg2) + 
#   geom_edge_link(aes(colour = lor), alpha = 1, show.legend = T) +
#   scale_edge_color_gradient2(name = "log(odds ratio)", low = "darkmagenta", mid = "white", high = "green1") +
#   geom_node_text(aes(label = species), size = 4, alpha = 0.5, repel = T) +
#   geom_node_point(aes(fill = community, size = degree), show.legend = F, shape = 21, alpha = 1) +
#   scale_colour_brewer(palette = "Set1", na.value = "grey50",) +
#   ggtitle('G2') +
#   theme_minimal() + th_no_axes( ) +
#   guides(fill = FALSE) +
#   guides(size = FALSE)
# plot(plt)

### CLUSTER LEADING EIGEN
# Remove imap clusters
# G1 <- G1 %>% select(-community)
# cle <- cluster_leading_eigen(G1)
# G1 <-(G1 %>% activate(nodes) 
#       %>% mutate(community = as.factor(membership(cle))))
# 
# lcle <- layout_with_fr(G1)
# G1 %>% activate(edges) %>% filter(lor > 0) %>% ggraph(layout = lcle) + 
#   geom_edge_link(aes(colour = lor), alpha = 1, show.legend = T) +
#   scale_edge_color_gradient2(name = "log(odds ratio)", low = "darkmagenta", mid = "white", high = "green1") +
#   geom_node_point(aes(fill = community, size = degree), show.legend = F, shape = 21, alpha = 1) +
#   # geom_node_text(aes(label=species), repel=TRUE, size =1.5,colour="black") +
#   ggtitle('G1 Cluster Leading Eigen') +
#   theme_minimal() + th_no_axes( ) +
#   guides(fill = FALSE) +
#   guides(size = FALSE)
# 
# cle_clusters <- as_tibble(G1 %>% activate(nodes))


