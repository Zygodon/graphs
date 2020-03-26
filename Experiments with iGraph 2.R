library(tidyverse)
library(igraph)
library(ggraph)
library("colorspace")

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

# Isolate between community and within community edges

df <- as_data_frame(G1a_sub, what = "both")
g <- graph_from_data_frame(df$edges, df$vertices, directed = F)
plot(g, vertex.label=NA, main = "g")
g_clust <- cluster_infomap(g)
plot(g_clust, g, vertex.label=NA, axes = F, rescale = T,
     ylim=c(-0.8, 0.8), xlim=c(-0.4, 0.4), asp = 1, vertex.label=NA, main = "g infomap")

clusters <- tibble(cluster = g_clust$membership, species = g_clust$names)
V(g)$community <- g_clust$membership
q4 <- qualitative_hcl(4, palette = "Dark 3")
V(g)$color <- q4[V(g)$community]

a1 <- as.data.frame(get.edgelist(g))

for (i in 1:length(E(g)))
{
  v1 <- a1$V1[i]
  v2 <- a1$V2[i]
  c1 <- inner_join(as_tibble(v1), clusters, by = c("value" = "species"))
  c2 <- inner_join(as_tibble(v2), clusters, by = c("value" = "species"))
  E(g)[i]$color <- ifelse(c1$cluster == c2$cluster, q4[c1$cluster], q4[4])
}
lo <- layout.fruchterman.reingold(g)
plot(g, vertex.label=NA, layout = lo, main = "test")

plt <- ggraph(g, layout=lo) + 
  geom_edge_link0(aes(color=E(g)$color), width=0.6) +
  geom_node_point(aes(color=V(g)$color), size=3)
plot(plt) # Note edge colours incorrect


