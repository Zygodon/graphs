
fp <- pairs %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
fc <- fp %>% group_by(species) %>% summarise(count = n())
pairs <- pairs %>% mutate(net_weights = (pairs$or - min(pairs$or))*2/max(pairs$or - min(pairs$or)))


net1 <- graph_from_data_frame(d = pairs, vertices = fc, directed = F)
l1 <- layout.fruchterman.reingold(net1)
# l <- layout_nicely(net1)
# l <- layout_with_lgl(net1)
#l <- layout_with_kk(net1)
# plot(net1, layout=l, vertex.size = 2, vertex.label=NA)
V(net1)$size <- 2 + V(net1)$count
E(net1)$weight <- pairs$net_weights
E(net1)$width <- abs(pairs$or*1.5)
E(net1)$color <- ifelse(pairs$or > 0, "green", "grey")
plot(net1, vertex.label=NA, layout = l1)

# largest_cliques(net1) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net1))
vcol[unlist(largest_cliques(net1))] <- "gold"
plot(as.undirected(net1), vertex.label=NA, vertex.color=vcol, layout = l1)


# cliques(net1) # list of cliques       
# sapply(cliques(net.sym), length) # clique sizes

largest_cliques(net1) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net1))
vcol[unlist(largest_cliques(net1))] <- "gold"
plot(as.undirected(net1), vertex.label="", vertex.color=vcol)
# plot(net1, vertex.color=vcol)

wc <- cluster_walktrap(net1, weights = E(net1)$weight)
modularity(wc)
plot(wc, net1, vertex.label=NA, layout = l)


# Friends
friends <- pairs %>% filter(or > 0.0)
frp <- friends %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
frc <- frp %>% group_by(species) %>% summarise(count = n())

net2 <- graph_from_data_frame(d = friends, vertices = fc, directed = F)
V(net2)$size <- 2 + V(net1)$count
E(net2)$weight <- friends$net_weights
E(net2)$width <- abs(friends$or*1.5)
E(net2)$color <- ifelse(friends$or > 0, "green", "grey")
# largest_cliques(net2) # cliques with max number of nodes
vcol <- rep("grey80", vcount(net1))
vcol[unlist(largest_cliques(net2))] <- "gold"
plot(net2, vertex.label=NA, vertex.color=vcol, layout = l1)

# sink("cliques.txt")
# largest_cliques(net2) # cliques with max number of nodes
# sink()

wc <- cluster_walktrap(net2, weights = E(net2)$weight)
plot(wc, net2, vertex.label = NA, layout = l1)

im <- cluster_infomap(net2)
plot(im, net2, vertex.label = NA, layout = l1)

sink("infomap.txt")
communities(im)
sink()

sink("degree_friends.txt")
degree(net2, mode="all")
sink()

deg <- degree(net2, mode="all")

