library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(DT)
library(blockmodels)

edges <- read_csv("edges.csv")
# Get a vertices list
vertices <- as_tibble(unique(unlist(c(edges$A, edges$B), " ")))

# G0: positive LOR only
G0 <- (tbl_graph(nodes = vertices, edges = edges, directed = F)
       %>% activate(edges)
       %>% filter(p_val < 0.05)
       # %>% filter(!(from == to))
       %>% filter(lor > 0)
       # %>% mutate(weight = lor)
       %>% activate(nodes)
       %>% mutate(degree = degree(., V(.)))
       %>% filter(degree > 5))

# Hairball preview
# G1 <-  as_tbl_graph(G0)
# Remove isolated nodes.
# G1 <- G1 %>% activate(nodes) %>% mutate(degree = degree(G1)) %>% filter(degree > 0)
G0 %>% activate(edges) %>% ggraph(layout = layout.fruchterman.reingold(G0)) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(shape = 21, alpha = 1) +
  ggtitle('G0 SBM Poisson preview') +
  theme_minimal()

# BLOCK MODEL
M <- as_adj(G0, type = "both", sparse = F)
my_model <- BM_bernoulli("SBM",M )
my_model$estimate()
model <- which.max(my_model$ICL) # Print which model max
# Get membership from the best model
mmZ <- my_model$memberships[[which.max(my_model$ICL)]]$Z
sbm_comms <- apply(mmZ, 1, which.max)
# Add node membership to the graph
G0 <-(G0 %>% activate(nodes) 
         %>% mutate(sbm_comm = as.factor(sbm_comms)))
# Add community membership to edges within each community
G0 <- (G0 %>% activate(edges)
       %>% mutate(sbm_comm = ifelse(.N()$sbm_comm[from] == .N()$sbm_comm[to], .N()$sbm_comm[from], NA)))

# Hairball with communities
lo <- layout_with_fr(G0)
G0 %>% activate(edges) %>% ggraph(layout = lo) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(aes(fill = sbm_comm, size = degree), show.legend = T, shape = 21, alpha = 1) +
  scale_fill_brewer(palette = "Dark2", na.value = "grey50") +
  ggtitle('G0 SBM Bernouilli, no weights, positive LOR') +
  theme_minimal() + th_no_axes( ) +
  guides(fill = guide_legend(override.aes = list(size=5))) # +

# Sort nodes by community size, SBM community 
# and species name (value)prior to looking at block matrix
no0 <- G0 %>% activate(nodes) %>% as_tibble(.)
community_size <- no0 %>% group_by(sbm_comm) %>% summarise(size = n())
G0 <- G0 %>% activate(nodes) %>% left_join(community_size, by = "sbm_comm")
G0 <- G0 %>% activate(nodes) %>% arrange(desc(size), sbm_comm, value)
# Re-do no0, sorted
no0 <- G0 %>% activate(nodes) %>% as_tibble(.)
# Convenient tibble
eg0 <- G0 %>% activate(edges) %>% as_tibble(.)

# Matrix
# This plot shows edges between communities as NA
pal <- brewer.pal(model, "Dark2")
no0 <- no0 %>% mutate(axis_colour = pal[sbm_comm])
# X, Y coordinates for text label reporting ICL
text_x <- length(row.names(no0)) - 10
text_y <-  -2
# Horizontal and vertical lines between communities
hvlines <- (no0 %>% group_by(sbm_comm)
            %>% summarise(n = n())
            %>% arrange(desc(n))
            %>% mutate(cs = cumsum(n) + 0.5))
                

ggraph(G0, 'matrix', sort.by = NULL) + 
  geom_edge_point(aes(colour = as.factor(sbm_comm)), mirror = TRUE) +
  guides(edge_colour = guide_legend(title = "community", override.aes = list(edge_size = 4))) +
  geom_vline(xintercept = c(0, hvlines$cs), alpha = 0.5, colour = "grey") +
  geom_hline(yintercept = c(0, hvlines$cs), alpha = 0.5, colour = "grey") +
  scale_edge_colour_brewer(palette = "Dark2", na.value = "grey50") +
  scale_y_reverse(breaks = seq(1, length(row.names(no0)), by = 1), labels = no0$value, "from") +
  scale_x_continuous(breaks = seq(1, length(row.names(no0)), by = 1), labels = no0$value, "to") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4.5, angle = 90, colour = no0$axis_colour, 
                                   face = 'bold', hjust = 1)) +
  theme(axis.text.y = element_text(size = 4.5, colour = no0$axis_colour, face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('G0 SBM Bernouilli, no weights, positive LOR') +
  annotate("text", label = paste("ICL = ", as.integer(max(my_model$ICL))), 
           x = text_x, y = text_y, size = 3, colour = "black")

# lggsave("SBM_Bernouilli_2.jpg", width = 20, height = 20, units = "cm")

  
