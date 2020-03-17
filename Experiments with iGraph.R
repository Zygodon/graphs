
# Libraries
library("RMySQL")
library(tidyverse)
library(igraph)
library(ggpmisc)
library(DT)
library(gridExtra)

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

# General SQL query
query <- function(q)
{
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

# Load the database
GetTheData <-  function()
{
  # GET DATA FROM DB
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  
  
  q <- sprintf('select assembly_id, assembly_name, quadrat_count, community, quadrat_id, quadrat_size, visit_date, records_id, species.species_id, 
    species.species_name from assemblies
      join quadrats on quadrats.assembly_id = assemblies_id
      join visit_dates on quadrats.vd_id = visit_dates.vds_id
      join records on records.quadrat_id = quadrats_id
      join species on species.species_id = records.species_id
      # Two assemblies have 0 quadrat count; exclude A.capillaris_stolonifera;
      # exclude some odd assemblies with no assigned community
    where quadrat_count = 5 and species.species_id != 4 and community is not null
    and quadrat_size = "2x2";') 
  # NOTE: this extract includes "MG5", i.e. some MG5 communities where 
  # the team have not decided
  # on a sub-group.
  
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

JointContingency <- function(d, A, B) # quadrat/assembly data, species_name, species_name
{
  d <- d %>% rename("id" = 1)
  A <- A[1]
  B <- B[1]
  As <- d %>% filter(species_name == A)
  Bs <- d %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "id")
  # Get all the assembly ids, including ones with neither A nor B
  q <- d %>% ungroup() %>% select(id) %>% distinct(id)
  j2 <- (left_join(q, j1, by = "id") 
         # NOTE: column length q >= column length j1
         %>% mutate(X1Y1 = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(X1Y0 = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(X0Y1 = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(X0Y0 = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  return(s)
}

# OddsRatio <-  function(s) # Joint contingency (colSums)
# {
#   # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
#   se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
#   o_r <- log((s[1]*s[4])/(s[2]*s[3])) # returning LOG OR
#   ci_low <- o_r - 1.96*se
#   ci_high <- o_r + 1.96*se
#   retval <- c(o_r, ci_low, ci_high)
#   names(retval) <- c("odds_ratio", "ci_low", "ci_high")
#   return(retval)
# }
# 
# Sverticeshi <-  function(jc)
# {
#   x <- matrix(unlist(jc), ncol = 2, nrow = 2, byrow = T)
#   ifelse(chisq.test(x)$p.value < 0.05, "yes", "no")
# }
######### END FUNCTIONS

the_data <- GetTheData()
assembly_count <- unlist(the_data %>% distinct(assembly_id) %>% summarise(n = n()))

# Make the edges 
assembly_data <- (the_data %>% select(assembly_id, species_name) 
                  %>% distinct()
                  %>% group_by(species_name, assembly_id))

edges <- (full_join(assembly_data, assembly_data, by = "assembly_id")
          %>% rename(A = species_name.x)
          %>% rename(B = species_name.y)
          %>% group_by(A, B) 
          %>% summarise(n = n())
          # %>% select(-n)
          %>% filter(A != B))
# # Get a vertices list
edge_pivot <- edges %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
vertices <- edge_pivot %>% group_by(species) %>% summarise(count = n()) %>% select(-count) #nodes; vertices
rm(edge_pivot)

G0 <- graph_from_data_frame(d = edges, vertices = vertices, directed = F)
# Remove loops and duplicate edges
G0 <- simplify(G0, edge.attr.comb = "first")
candidate_edges <- as_tibble(as_data_frame(G0, what="edges")) %>% rename(A = from, B = to)
rm(edges, G0) # G0 not suitable for analysis ...

# For each species pair, get the odds ratio.
aor <- tibble(
  A = candidate_edges$A,
  B = candidate_edges$B,
  jc1 = 0,
  jc2 = 0,
  jc3 = 0,
  jc4 = 0,
  or = 0,
  pval = 0)

for (i in seq_along(row.names(candidate_edges)))
{
  s <- JointContingency(assembly_data, candidate_edges$A[i], candidate_edges$B[i])
  aor[i, 3:6] <- s[1:4]
}
# Remove 0 entries in the contingency tables
aor <-(aor %>% filter(jc1 != 0)
       %>% filter(jc2 != 0)
       %>% filter(jc3 != 0)
       %>% filter(jc4 != 0))

for (i in seq_along(row.names(aor)))
{
  s <- aor[i, 3:6]
  ft <- fisher.test(matrix(s, nrow = 2, ncol = 2, byrow = T) ,alternative = "greater")
  aor[i, 7] <- ft$estimate  # or
  aor[i, 8] <- ft$p.value
}
# 
# # aor <- aor %>% mutate(or = ifelse(is.infinite(or), 300, or)) %>% mutate(lor = log(or))
# 
# aor <- aor %>% mutate(lor = log(or))
# # aor <- (aor %>% filter(!is.na(or)) # Actually there are no NAs with fisher estimate
#             # remove cases, A with no B, B with no A. 1770 of these
# #            %>% filter(!is.infinite(or))) # infinites come if A has no B & vice versa
# 
# # Scale 0 =< lor < 2 for net weights (used as a scale factor, so most -ve lor is weight 0
# # so should not be represented as an edge)
# # aor <- aor %>% mutate(net_weights = (aor$lor - min(aor$lor))*2/max(aor$lor - min(aor$lor)))
# 
# # Check on the lor histogram
# fig2a <- ggplot(aor, aes(lor))  + 
#   geom_histogram(aes(y = ..density..), binwidth = 0.25, colour = "black") + 
#   stat_function(fun = dnorm, args = list(mean = mean(aor$lor), sd = sd(aor$lor)), colour = "green") +
#   geom_vline(xintercept = 0, colour = "red") +
#   xlim(-3, 6) +
#   ylim(0.0, 0.4) +
#   labs(title ="Figure 2a", x = "log(odds ratio)") 
# plot(fig2a)
# 
# # Get a vertices list from the species names still in aor
# fp <- aor %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
# vertices <- fp %>% group_by(species) %>% summarise(count = n()) #nodes; vertices
# rm(fp)
# ## EDGE TRIMMING BY "VALUE"
# 
# aor <- aor %>% mutate(lor2 = lor^2)
# # Ordered sequence of values to control edge removal
# values <- aor %>%  select(lor2) %>% arrange(lor2) %>% filter(lor2 > 0.0) %>% distinct(lor2)
# 
# net1 <- graph_from_data_frame(d = aor, vertices = vertices, directed = F)
# # l1 <- layout.fruchterman.reingold(net1)
# # V(net1)$size <- 10 * V(net1)$count/max(V(net1)$count)
# V(net1)$size <- 1
# # E(net1)$weight <- aor$net_weights
# # E(net1)$width <- 0.5 * abs(aor$lor)
# E(net1)$color <- ifelse(aor$lor > 0, "green", "grey")
# #plot(net1, vertex.label=NA, layout = l1)
# plot(net1, vertex.label=NA)
# 
# 
# mods <- tibble(mod = rep(0.0, length(row.names(values))),
#                edges = 0, 
#                vertices = 0, 
#                i = 0)
# 
# for (i in seq_along(row.names(values)))
# {
#   net2 <- delete_edges(net1, which(E(net1)$lor2 < values$lor2[i]))
#   isolated <-  which(degree(net2)==0)
#   net2 <-  delete.vertices(net2, isolated)
#   mods$mod[i] <- modularity(cluster_infomap(net2))
#   mods$edges[i] <- length(E(net2))
#   mods$vertices[i] <- length(V(net2))
#   mods$i[i] <- i
# }
# 
# ggplot(mods, aes(i, mod)) +
#   geom_point()
# 
# net_low <- delete_edges(net1, E(net1)[which(E(net1)$lor2 < values$lor2[3000])])
# im_low <- cluster_infomap(net_low)
# modularity(im_low)
# plot(im_low, net_low, vertex.label=NA)
# 
# net_high <- delete_edges(net1, E(net1)[which(E(net1)$lor2 < values$lor2[3200])])
# im_high <- cluster_infomap(net_high)
# modularity(im_high)
# plot(im_high, net_high, vertex.label=NA)
# 
# sink("communities im_low  2020-03-13.txt")
# print(communities(im_low))
# sink()



### EDGE TRIMMING BY lor QUANTILES
# net1 <- graph_from_data_frame(d = aor, vertices = vertices, directed = F)
# l1 <- layout.fruchterman.reingold(net1)
# # V(net1)$size <- 10 * V(net1)$count/max(V(net1)$count)
# V(net1)$size <- 1
# E(net1)$weight <- aor$net_weights
# # E(net1)$width <- 0.5 * abs(aor$lor)
# # E(net1)$color <- ifelse(aor$lor > 0, "green", "grey")
# plot(net1, vertex.label=NA, layout = l1)
# 
# 
# im <- cluster_infomap(net1)
# modularity(im)
# plot(im, net1, vertex.label = NA, layout = l1)

# Remove edges by quantile, starting by removing all negative associations

# q <- quantile(E(net1)$weight, probs = seq(0.5, 1, 0.01))
# e_count <- length(E(net1))
# mods <- tibble(mod = rep(0.0, length(q)),
#                edges = 0, 
#                vertices = 0, 
#                i = 0)
# 
# for (i in seq_along(row.names(mods)))
# {
#   net2 <- delete_edges(net1, E(net1)[which(E(net1)$weight < q[i])])
#   isolated <-  which(degree(net2)==0)
#   net2 <-  delete.vertices(net2, isolated)
#   mods$mod[i] <- modularity(cluster_infomap(net2))
#   mods$edges[i] <- length(E(net2))
#   mods$vertices[i] <- length(V(net2))
#   mods$i[i] <- i
# }
# 
# ggplot(mods, aes(i, mod)) +
#   geom_point()
# 
# ggplot(mods, aes(i, vertices)) +
#   geom_point()
# 
# 
# net3 <- delete_edges(net1, E(net1)[which(E(net1)$weight < q[45])])
# im <- cluster_infomap(net3)
# modularity(im)


# ### Vertex trimming: Does not change modularity (stays at 0)
# 
# vsA <- (aor %>% select(A, B, lor) 
#           %>% mutate(lor2 = lor^2)
#           %>% group_by(A) 
#           %>% summarise(s2A = sum(lor2)))
# vsB <- (aor %>% select(A, B, lor) 
#         %>% mutate(lor2 = lor^2)
#         %>% group_by(B) 
#         %>% summarise(s2B = sum(lor2)))
# 
# t <- full_join(vsA, vsB, by = c("A" = "B"))
# t <- (t %>% replace_na(list(s2A = 0, s2B = 0))
#       %>% mutate(lor2 = s2A + s2B)
#       %>% select(species = A, value = lor2))
# 
# values <- t %>%  arrange(value)
# v_mods <- tibble(mod = rep(0.0, length(values$value)),
#                edges = 0, 
#                vertices = 0, 
#                i = 0)
# 
# g1 <- graph_from_data_frame(d = aor, vertices = t, directed = F)
# 
# for (i in seq_along(row.names(t)))
# {
#   g2 <- delete_vertices(g1, which(V(g1)$value < values$value[i]))
#   isolated <-  which(degree(g2)==0)
#   g2 <-  delete.vertices(g2, isolated)
#   v_mods$mod[i] <- modularity(cluster_infomap(g2))
#   v_mods$edges[i] <- length(E(g2))
#   v_mods$vertices[i] <- length(V(g2))
#   v_mods$i[i] <- i
# }
# 
# ggplot(v_mods, aes(i, mod)) +
#   geom_point()
# 
# ggplot(v_mods, aes(i, vertices)) +
#   geom_point()


######### SOME NETWORK METRICS ORIGINALLY MADE ON G0
# # Let's have a look at G0 before going any farther.
# plot(G0, vertex.label=NA)
# # And get some data about it
# degG0 <- tibble(degree(G0))
# vertices <- bind_cols(vertices, degG0) 
# vertices <- vertices %>% rename(degree = "degree(G0)")
# rm(degG0)
# 
# # Check on vertex degree histogram
# p1 <- ggplot(vertices, aes(degree))  +
#   geom_histogram() +
#   labs(title ="Node degrees", x = "degree")
# plot(p1)
# 
# hub_score_G0 <- tibble(hub_score(G0)$vector)
# # The hub scores of the vertices are defined as the principal eigenvector of AAT,
# # where A is the adjacency matrix of the graph.
# vertices <- bind_cols(vertices, hub_score_G0)
# vertices <- vertices %>% rename(hub_score = `hub_score(G0)$vector`)
# 
# # vertex_connectivity: the minimum number of vertices to remove to make the graph
# # "not strongly connected".
# vcG0 <- vertex_connectivity(G0)# =16
# # Transitivity: "Transitivity measures the probability that 
# # the adjacent vertices of a vertex are connected."
# trans_G0 <- transitivity(G0, type = "global") # 0.632

# del <- vertices %>% filter(hub_score > 0.7)
# G1 <- delete_vertices(G0, del$species)
# plot(G1, vertex.label=NA)
# imG1 <- cluster_infomap(G1)
# modularity(imG1)
# transitivity(G1)
# plot(imG1, G1, vertex.label=NA, main = "hub score threshold 0.7")



