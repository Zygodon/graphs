# STATISTICAL BLOCK MODEL USING POISSON DISTRIBUTION 
# ON GRAPH EDGES - POSITIVE AND NEGATIVE LOG(ODDS RATIO)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(blockmodels)

# FUNCTIONS

posxy <-  function(hvlines, comms_order)
{
  k <- length(hvlines)
  locs <- c(0, hvlines)
  # ys <- matrix(rep(0, 4*k), 2*k, 2)         # Always 4 corners to a rectangle
  ys <- matrix(rep(0, 2*length(locs)), length(locs), 2)         # Always 4 corners to a rectangle
  # for(i in 1:4) 
  for(i in 1:length(locs)) 
  {
    ys[i,] <- rep(locs[i], 2) # y coords top and bottom hvlines
  }
  y2 <- matrix(rep(0, 4 * k), k, 4)       # matrix, zeros
  # Fill y2
  for (i in 1:k)
  {
    y2[i,] <- c(ys[i,], ys[i + 1,]) # y coords for each group
  }
  x2 <- y2[, c(1,4,3,2)]                  # Swap columns 2, 4
  y3 <- as.vector(t(y2))
  x3 <- as.vector(t(x2))
  
  pos <- data.frame(
    comm = rep(comms_order, each = 4),
    x = x3,
    y = y3)
  return(pos)
}

# MAIN

# Edges stored in file. See Stand occupancy.R for code
# Key: A, species-from; B, species-to; a,b,c,d, entries in the contingency tabel
# p_val, significance level, Fisher's Exact Text for association;
# or, odds ratio
# lor, logarithm of the odds ratio
edges <- read_csv("edges.csv")
# Get a vertices list
vertices <- as_tibble(unique(unlist(c(edges$A, edges$B), " ")))

mu <- mean(edges$lor)
stdev <- sd(edges$lor)
# # Plot parameters, show in block matrix plot text
p_lim <- 0.05
text_size = 5

# G0: the basic graph
edges1 <- (edges %>% mutate(sp_from = A, sp_to = B)
           %>% filter(lor < 0)
           %>% filter(p_val < p_lim)
           %>% filter((lor - (1.95*ase)>= 0) | (lor + (1.95*ase))<= 0))

vertices1 <- edges1 %>% pivot_longer(cols = c(A, B), values_to = "species") %>% distinct(species)
G0 <- (tbl_graph(nodes = vertices1, edges = edges1, directed = F)
              %>% mutate(deg = centrality_degree()))
ggplot(data = enframe(igraph::degree(G0)), aes(x = value)) + 
  geom_density(size = 1)
G1 <- G0 %>% filter(deg >= 7)

G2 <- (G0 %>% activate(nodes) 
          %>%  mutate(edge_node = node_is_adjacent(V(G0)[which(V(G0)$species %in% V(G1)$species)]))
          %>% filter((edge_node == T)))

# Hairball preview
G2 %>% activate(edges) %>% ggraph(layout = layout.fruchterman.reingold(G2)) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(shape = 21, alpha = 1) +
  ggtitle('G2 preview') +
  theme_minimal()

# BLOCK MODEL
M <- as_adj(G0, type = "both", sparse = F)
C <- as_adj(G0, type = "both", attr = "lor", sparse = F)
D <- as_adj(G2, type = "both", attr = "a", sparse = F)

# my_model <- BM_bernoulli("SBM_sym",M )
# my_model <- BM_bernoulli_covariates_fast("SBM",M,C)
# my_model <- my_model <- BM_gaussian("SBM",C )
my_model <- BM_poisson("SBM_sym",D)
# my_model <- BM_poisson_covariates("SBM_sym", D, C)
# my_model <- BM_poisson_covariates("SBM", D, C)

my_model$estimate()
# NOTE: ICL Integrated Completed Likelihood: figure of merit for the multiple
# models explored during estimation
model <- 6 # which.max(my_model$ICL) # Select the model to use
# Print C matrix
print(my_model$model_parameters[model])
# Get membership from the best model
# mmZ <- my_model$memberships[[which.max(my_model$ICL)]]$Z
mmZ <- my_model$memberships[[model]]$Z
sbm_comms <- apply(mmZ, 1, which.max)

# Add node membership to the graph
G2 <-(G2 %>% activate(nodes) 
      %>% mutate(sbm_comm = as.factor(sbm_comms)))
# Add community membership to edges within each community
# NOTE: edges BETWEEN communities flagged as NA
G2 <- (G2 %>% activate(edges)
       %>% mutate(sbm_comm = ifelse(.N()$sbm_comm[from] == .N()$sbm_comm[to], .N()$sbm_comm[from], NA)))

# Hairball with communities
lo <- layout_with_fr(G2)
G2 %>% activate(edges) %>% ggraph(layout = lo) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(aes(fill = sbm_comm), size = 4, show.legend = T, shape = 21, alpha = 1) +
  scale_fill_brewer(palette = "Dark2", na.value = "grey50") +
  ggtitle('SBM Poisson') +
  theme_minimal() + th_no_axes( ) +
  guides(fill = guide_legend(override.aes = list(size=5))) 
# ggsave("SBM_Poisson_hairball.jpg", width = 10, height = 10, units = "cm")

# Sort nodes by community size, SBM community and
# species name (value) prior to looking at block matrix
no2 <- G2 %>% activate(nodes) %>% as_tibble(.)
community_size <- no2 %>% group_by(sbm_comm) %>% summarise(size = n())
G2 <- G2 %>% activate(nodes) %>% left_join(community_size, by = "sbm_comm")
G2 <- G2 %>% activate(nodes) %>% arrange(desc(size), sbm_comm, species)
# Re-do no0, sorted
no2 <- G2 %>% activate(nodes) %>% as_tibble(.)
# Convenient tibble
eg2 <- G2 %>% activate(edges) %>% as_tibble(.)

# Matrix
# This plot shows edges between communities as NA
pal <- brewer.pal(model, "Dark2") # model: the selected SBM
no2 <- no2 %>% mutate(axis_colour = pal[sbm_comm])
# X, Y coordinates for text label reporting ICL
text_x <- length(row.names(no2)) - 20
text_y <-  -3
# Text to show parameters on plot
# txt <- paste(paste("ICL =", as.integer(max(my_model$ICL)), sep = " "), 
#              paste("p_val <", p_lim, sep = " "),
#              paste("LOR low =", format(lor_low, digits = 2),
#                    "LOR high =", format(lor_high, digits = 2), sep = " "), sep = "\n")
# Horizontal and vertical lines between communities
hvlines <- (no2 %>% group_by(sbm_comm)
            %>% summarise(n = n())
            %>% arrange(desc(n))
            %>% mutate(cs = cumsum(n) + 0.5))
# Rectangles
rectangles <- posxy(hvlines$cs, hvlines$sbm_comm)
pal <- brewer.pal(7, "Dark2")
pal1 <- pal[as.numeric(levels(as.factor(eg2$sbm_comm)))]

ggraph(G2, 'matrix', sort.by = NULL) + 
  geom_polygon(data = rectangles, aes(x = x, y = y, fill = as.factor(comm), group = comm), alpha = 1) +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(title = "community", override.aes = list(alpha = 1))) +
  geom_edge_point(aes(colour = lor, size = a), edge_alpha = 0.6, mirror = TRUE) +
  scale_edge_colour_gradient2(low = "black",
                              mid = "#f7f7f7", high = "#af8dc3", midpoint = 2) +
  geom_edge_point(edge_shape = 3, edge_size = 0.1, edge_alpha = 1, mirror = TRUE) +
  guides(edge_size = guide_legend(title = "instances")) +
  geom_vline(xintercept = c(0, hvlines$cs), alpha = 0.5, colour = "grey") +
  geom_hline(yintercept = c(0, hvlines$cs), alpha = 0.5, colour = "grey") +
  ## scale_edge_alpha(trans = 'reverse') + # Reverse the alpha scale for negative LOR
  scale_y_reverse(breaks = seq(1, length(row.names(no2)), by = 1), labels = no2$species, "from") +
  scale_x_continuous(breaks = seq(1, length(row.names(no2)), by = 1), labels = no2$species, "to") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(size = text_size, angle = 90, colour = no2$axis_colour, 
                                   face = 'bold', hjust = 1)) +
  theme(axis.text.y = element_text(size = text_size, colour = no2$axis_colour, face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.box = "vertical") +
  ggtitle('SBM Poisson') # +
 # annotate("text",
 #    label = txt,
 #    x = text_x,
 #    y = text_y,
 #    size = 2,
 #    colour = "black")

# ggsave("SBM_Poisson__LOR.jpg", width = 20, height = 20, units = "cm")
