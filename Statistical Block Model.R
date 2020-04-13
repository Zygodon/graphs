library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(DT)
library(blockmodels)


# site_occ <- read_csv("stand occupancy.csv")
edges <- read_csv("edges.csv")
# Get a vertices list
edge_pivot <- edges %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
vertices <- (edge_pivot %>% group_by(species) 
             %>% summarise(count = n()) 
             %>% select(-count) #nodes; vertices
             %>% rename(name = species))
rm(edge_pivot)

# G0: positive LOR only
G0 <- (tbl_graph(nodes = vertices, edges = edges, directed = F)
       %>% activate(edges)
       %>% filter(p_val < 0.05)
       %>% filter(!(from == to))
       %>% filter(lor > 0)
       %>% mutate(weight = lor)
       %>% activate(nodes)
       %>% mutate(degree = degree(., V(.)))
       %>% filter(degree > 5))

# BLOCK MODEL
M <- as_adj(G0, type = "both", sparse = F)
my_model <- BM_bernoulli("SBM",M )
my_model$estimate()
which.max(my_model$ICL)
mmZ <- my_model$memberships[[which.max(my_model$ICL)]]$Z
sbm_comms <- apply(mmZ, 1, which.max)

G0 <-(G0 %>% activate(nodes) 
      %>% mutate(sbm_comm = as.factor(sbm_comms)))
# Identify edges WITHIN communities
G0 <- (G0 %>% activate(edges)
       %>% mutate(sbm_comm = ifelse(.N()$sbm_comm[from] == .N()$sbm_comm[to], .N()$sbm_comm[from], NA)))
eg0 <- G0 %>% activate(edges) %>% as_tibble(.)
no0 <- G0 %>% activate(nodes) %>% as_tibble(.)

ggraph(G0, 'matrix', sort.by = node_rank_leafsort()) + 
  # ggraph(G0, 'matrix', sort.by = node_rank_spectral()) + 
  geom_edge_point(aes(colour = as.factor(sbm_comm)), mirror = TRUE) + 
  theme(legend.position = 'bottom')

# Sort nodes by SBM community prior to looking at block matrix
G1 <- G0 %>% activate(nodes) %>% arrange(sbm_comm)

# Hairball
sp <- V(G1)$name
lo <- layout_with_fr(G1)
G1 %>% activate(edges) %>% ggraph(layout = lo) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(aes(fill = sbm_comm, size = degree), show.legend = T, shape = 21, alpha = 1) +
  scale_fill_brewer(palette = "Dark2", na.value = "grey50") +
  ggtitle('G0 SBM Bernouilli, no weights, positive LOR') +
  theme_minimal() + th_no_axes( ) +
  guides(fill = guide_legend(override.aes = list(size=5))) # +

# Matrix
# This plot shows edges between communities as NA
pal <- brewer.pal(6, "Dark2")
no1 <- G1 %>% activate(nodes) %>% as_tibble(.)
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
  scale_y_reverse(breaks = seq(1, 79, by = 1), labels = no1$name, "from") +
  scale_x_continuous(breaks = seq(1, 79, by = 1), labels = no1$name, "to") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4.5, angle = 90, colour = no1$axis_colour, 
                                   face = 'bold', hjust = 1)) +
  theme(axis.text.y = element_text(size = 4.5, colour = no1$axis_colour, face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('G1 SBM Bernouilli, no weights, positive LOR') +
  annotate("text", label = paste("ICL = ", as.integer(max(my_model$ICL))), 
           x = text_x, y = text_y, size = 3, colour = "black")

ggsave("SBM_1.jpg", width = 20, height = 20, units = "cm")

  



# cle <- cluster_leading_eigen(G0, weights = NA)
# G0 <-(G0 %>% activate(nodes) 
#       %>% mutate(community = as.factor(membership(cle))))
# cle_clusters <- as_tibble(G0 %>% activate(nodes))



# no0 <- no0 %>% mutate(sbm_comm = apply(mmZ, 1, which.max))





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
# lo <- layout_with_fr(G1)
# G1 %>% activate(edges) %>% filter(lor > 0) %>% ggraph(layout = lo) + 
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


