# STATISTICAL BLOCK MODEL USING BERNOUILLI/BINOMIAL DISTRIBUTION 
# ON GRAPH EDGES.
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

# Plot parameters, show in block matrix plot text
p_lim <- 0.005
deg <- 2

# G0: the basic graph
G0 <- (tbl_graph(nodes = vertices, edges = edges, directed = F)
       %>% activate(edges)
       %>% filter(p_val < p_lim)
       %>% filter(lor > 0)
       %>% activate(nodes)
       %>% mutate(degree = degree(., V(.)))
       %>% filter(degree > deg))
    
# Hairball preview
G0 %>% activate(edges) %>% ggraph(layout = layout.fruchterman.reingold(G0)) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(shape = 21, alpha = 1) +
  ggtitle('G0 preview') +
  theme_minimal()

# BLOCK MODEL
M <- as_adj(G0, type = "both", sparse = F)
my_model <- BM_bernoulli("SBM",M )
my_model$estimate()
# NOTE: ICL Integrated Completed Likelihood: figure of merit for the multiple
# models explored during estimation
model <- which.max(my_model$ICL) # Select the model to use
# Get membership from the best model
mmZ <- my_model$memberships[[which.max(my_model$ICL)]]$Z
sbm_comms <- apply(mmZ, 1, which.max)
# Add node membership to the graph
G0 <-(G0 %>% activate(nodes) 
      %>% mutate(sbm_comm = as.factor(sbm_comms)))
# Add community membership to edges within each community
# NOTE: edges BETWEEN communities flagged as NA
G0 <- (G0 %>% activate(edges)
       %>% mutate(sbm_comm = ifelse(.N()$sbm_comm[from] == .N()$sbm_comm[to], .N()$sbm_comm[from], NA)))

# Hairball with communities
lo <- layout_with_fr(G0)
G0 %>% activate(edges) %>% ggraph(layout = lo) + 
  geom_edge_link(colour = "black", alpha = 0.2) +
  geom_node_point(aes(fill = sbm_comm, size = degree), show.legend = T, shape = 21, alpha = 1) +
  scale_fill_brewer(palette = "Dark2", na.value = "grey50") +
  ggtitle('SBM Bernouilli') +
  theme_minimal() + th_no_axes( ) +
  guides(fill = guide_legend(override.aes = list(size=5))) 
# ggsave("SBM_Bernouilli_4_hairball.jpg", width = 10, height = 10, units = "cm")

# Sort nodes by community size, SBM community and
# species name (value) prior to looking at block matrix
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
pal <- brewer.pal(model, "Dark2") # model: the selected SBM
no0 <- no0 %>% mutate(axis_colour = pal[sbm_comm])
# X, Y coordinates for text label reporting ICL
text_x <- length(row.names(no0)) - 10
text_y <-  -2
# Text to show parameters on plot
txt <- paste(paste("ICL =", as.integer(max(my_model$ICL)), sep = " "), 
             paste("p_val <", p_lim, sep = " "),
             paste("degree >", deg, sep = " "), sep = "\n")
# Horizontal and vertical lines between communities
hvlines <- (no0 %>% group_by(sbm_comm)
            %>% summarise(n = n())
            %>% arrange(desc(n))
            %>% mutate(cs = cumsum(n) + 0.5))
# Rectangles
rectangles <- posxy(hvlines$cs, hvlines$sbm_comm)

ggraph(G0, 'matrix', sort.by = NULL) + 
  geom_polygon(data = rectangles, aes(x = x, y = y, fill = comm, group = comm), alpha = 0.2) +
  scale_fill_brewer(palette = "Dark2", na.value = "grey50") +
  guides(fill = guide_legend(title = "community", override.aes = list(alpha = 1))) +
  geom_edge_point(aes(colour = as.factor(sbm_comm), size = a, alpha = lor), mirror = TRUE) +
  guides(colour = FALSE) + #No guide for edge point colour
  geom_edge_point(edge_shape = 3, edge_size = 0.1, edge_alpha = 0.5, mirror = TRUE) +
  guides(edge_alpha = guide_legend(title = "log(odds ratio)", override.aes = list(edge_size = 4))) +
  guides(edge_size = guide_legend(title = "instances")) +
  geom_vline(xintercept = c(0, hvlines$cs), alpha = 0.5, colour = "grey") +
  geom_hline(yintercept = c(0, hvlines$cs), alpha = 0.5, colour = "grey") +
  scale_edge_colour_brewer(palette = "Dark2", na.value = "grey50", guide = F) +
  # scale_edge_alpha(trans = 'reverse') + # Reverse the alpha scale for negative LOR
  scale_y_reverse(breaks = seq(1, length(row.names(no0)), by = 1), labels = no0$value, "from") +
  scale_x_continuous(breaks = seq(1, length(row.names(no0)), by = 1), labels = no0$value, "to") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, colour = no0$axis_colour, 
                                   face = 'bold', hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = no0$axis_colour, face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.box = "vertical") +
  ggtitle('SBM Bernouilli positive LOR') +
  annotate("text", 
           label = txt, 
           x = text_x, 
           y = text_y, 
           size = 2, 
           colour = "black")

# ggsave("SBM_Bernouilli_4.jpg", width = 20, height = 20, units = "cm")




