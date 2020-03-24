library(tidyverse)
library(igraph)

site_occ <- read_csv("stand occupancy.csv")
edges <- read_csv("edges.csv")
# Get a vertices list
edge_pivot <- edges %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
vertices <- edge_pivot %>% group_by(species) %>% summarise(count = n()) %>% select(-count) #nodes; vertices
rm(edge_pivot)

# Check on the lor histogram
plt1 <- ggplot(edges, aes(lor))  +
  geom_histogram(aes(y = ..density..), binwidth = 0.25, colour = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(edges$lor), sd = sd(edges$lor)), colour = "green") +
  geom_vline(xintercept = 0, colour = "red") +
  xlim(-3, 6) +
  ylim(0.0, 0.4) +
  labs(title ="", x = "log(odds ratio)")
plot(plt1)
t.test(edges$lor, mu = 0) # p-value < 2.2e-16

G0 <- graph_from_data_frame(d = edges, vertices = vertices, directed = F)
l0 <- layout.fruchterman.reingold(G0, weights = edges$or)
E(G0)$weight <- edges$or
plot(G0, layout = l0, vertex.label=NA, main = "G0")


e1 <- edges %>% filter(p_val < 0.05)
G1 <- graph_from_data_frame(d = e1, vertices = vertices, directed = F)
E(G1)$color <- ifelse(e1$lor > 0, "green", "grey50")
E(G1)$weight <- NA #e1$lor
V(G1)$size = 1
V(G1)$color <- "tomato2"
isolated <-  which(degree(G1)==0)
G1 <-  delete.vertices(G1, isolated)
# l1 <- layout.fruchterman.reingold(G1) #, weights = e1$lor)
plot(G1, vertex.label=NA, main = "G1")
# plot(G1, layout = l1, vertex.label=NA)

# Concentrate on the positive associations:
e1a <-  e1 %>% filter(lor > 0) # Positive associations
G1a <- graph_from_data_frame(d = e1a, vertices = vertices, directed = F)
l1a <- layout.fruchterman.reingold(G1a, weights = e1a$lor)
isolated <-  which(degree(G1a)==0)
G1a <-  delete.vertices(G1a, isolated)
E(G1a)$weight <- e1a$or
plot(G1a, vertex.label=NA, main = "G1a")
h <- hub_score(G1a)$vector
h_score <- tibble(species = names(h), score = h)
inc_edges <- incident(G1a,  V(G1a)["Lathyrus_montanus"], mode="all")

# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(G1a))
ecol[inc_edges] <- "orange"
vcol <- rep("grey40", vcount(G1a))
vcol[V(G1a)$name=="Lathyrus_montanus"] <- "gold"
V(G1a)$size = degree(G1a)/5
plot(G1a, vertex.color=vcol, edge.color=ecol, vertex.label=NA, axes = F, rescale = T,
     ylim=c(-0.8, 0.8), xlim=c(-0.4, 0.4), asp = 1, vertex.label=NA, main = "G1a L_montanus")

im <- cluster_infomap(G1a)
plot(im, G1a, vertex.label=NA, axes = F, rescale = T,
     ylim=c(-0.8, 0.8), xlim=c(-0.4, 0.4), asp = 1, vertex.label=NA, main = "G1a infomap")

modularity(im)

# Keep communities with 10 or more members
# https://stackoverflow.com/questions/51271730/how-to-remove-small-communities-using-igraph-in-r
# hermidalc April 2019
im_keep_ids <- as.numeric(names(sizes(im)[sizes(im) >= 10]))
im_keep_v_idxs <- which(im$membership %in% im_keep_ids)

G1a_sub <- induced_subgraph(G1a, V(G1a)[im_keep_v_idxs])
# igraph has no direct functionality to subset community objects so hack it
im_sub <- im
im_sub$names <- im$names[im_keep_v_idxs]
im_sub$membership <- im$membership[im_keep_v_idxs]
im_sub$vcount <- length(im_sub$names)
im_sub$modularity <- modularity(G1a_sub, im_sub$membership, E(G1a_sub)$weight)
V(G1a_sub)$size = degree(G1a_sub)/2
plot(im_sub, G1a_sub, vertex.label=NA, axes = F, rescale = T,
     ylim=c(-1, 1), xlim=c(-0.4, 0.4), asp = 1, vertex.label=NA, main = "G1a_sub")
modularity(im_sub)

# This gives 4 communities: attempt SBM
e1a_sub <- as_tibble(as_data_frame(G1a_sub, what = c("edges")))
v1a_sub <- as_tibble(as_data_frame(G1a_sub, what = c("vertices")))
m <- tibble(membership(im_sub))
v1a_sub <- bind_cols(v1a_sub, m) %>% rename(group = 3)
rm(m)
# Note v1a_sub vertex size = degree/2
