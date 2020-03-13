
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
  # Get all the assemblies/quadrats, including ones with neither A nor B
  q <- d %>% distinct(id) 
  j2 <- (left_join(q, j1, by = "id") 
         # NOTE: column length q >= column length j1
         %>% mutate(X1Y1 = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(X1Y0 = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(X0Y1 = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(X0Y0 = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  return(s)
}

OddsRatio <-  function(s) # Joint contingency (colSums)
{
  # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  o_r <- log((s[1]*s[4])/(s[2]*s[3])) # returning LOG OR
  ci_low <- o_r - 1.96*se
  ci_high <- o_r + 1.96*se
  retval <- c(o_r, ci_low, ci_high)
  names(retval) <- c("odds_ratio", "ci_low", "ci_high")
  return(retval)
}

SfChi <-  function(jc)
{
  x <- matrix(unlist(jc), ncol = 2, nrow = 2, byrow = T)
  ifelse(chisq.test(x)$p.value < 0.05, "yes", "no")
}
######### END FUNCTIONS

the_data <- GetTheData()
# Restrict analysis to species with more than 20 hits
sp_count <- (the_data 
             %>% group_by(species_id) 
             %>% summarise(count = n()))
n_sp <- length(row.names(sp_count))
n_pairs <- n_sp * (n_sp-1)
assembly_count <- the_data %>% distinct(assembly_id) %>% summarise(n = n())

# # DON'T REDUCE THE DATA!

# Make the edges
assembly_data <- the_data %>% select(assembly_id, species_name) %>% distinct()
quadrat_data <- the_data %>% select(quadrat_id, species_name)
edges <- (full_join(quadrat_data, quadrat_data, by = "quadrat_id")
          %>% rename(from = species_name.x)
          %>% rename(to = species_name.y)
          %>%  select(-quadrat_id)
          %>% group_by(from, to) 
          %>% summarise(share_2x2 = n())
          %>% filter(from != to))
nodes <- distinct(edges, from)
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
n2 <- simplify(net, edge.attr.comb = "first")
edges <- as_data_frame(n2, what="edges")
rm(net, n2)

# For each species pair, get the assembly odds ratio.
aor <- tibble(
  A = edges$from,
  B = edges$to,
  jc1 = 0,
  jc2 = 0,
  jc3 = 0,
  jc4 = 0,
  or = 0)

for (i in seq_along(row.names(edges)))
{
  s <- JointContingency(assembly_data, edges$from[i], edges$to[i])
  aor[i, 3:6] <- s[1:4]
  ft <- fisher.test(matrix(s, nrow = 2, ncol = 2, byrow = T) ,alternative = "greater")
  aor[i, 7] <- ft$estimate
}

aor <- aor %>% mutate(lor = log(or))
aor <- (aor %>% filter(!is.na(or))
        %>% filter(!is.infinite(or))
        %>% filter(!is.infinite(lor)))

fig2a <- ggplot(aor, aes(lor))  + 
  geom_histogram(aes(y = ..density..), binwidth = 0.25, colour = "black") + 
  stat_function(fun = dnorm, args = list(mean = mean(aor$lor), sd = sd(aor$lor)), colour = "green") +
  geom_vline(xintercept = 0, colour = "red") +
  xlim(-3, 6) +
  ylim(0.0, 0.4) +
  labs(title ="Figure 2a", x = "log(odds ratio)") 
plot(fig2a)

aor <- aor %>% mutate(net_weights = (aor$lor - min(aor$lor))*2/max(aor$lor - min(aor$lor)))

fp <- aor %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
fc <- fp %>% group_by(species) %>% summarise(count = n()) #nodes; vertices

net1 <- graph_from_data_frame(d = aor, vertices = fc, directed = F)
l1 <- layout.fruchterman.reingold(net1)
# V(net1)$size <- 10 * V(net1)$count/max(V(net1)$count)
V(net1)$size <- 1
E(net1)$weight <- aor$net_weights
# E(net1)$width <- 0.5 * abs(aor$lor)
# E(net1)$color <- ifelse(aor$lor > 0, "green", "grey")
plot(net1, vertex.label=NA, layout = l1)


im <- cluster_infomap(net1)
modularity(im)
plot(im, net1, vertex.label = NA, layout = l1)

# Remove edges by quantile, starting by removing all negative associations

q <- quantile(E(net1)$weight, probs = seq(0.5, 1, 0.01))
e_count <- length(E(net1))
mods <- tibble(mod = rep(0.0, length(q)),
               edges = 0, 
               vertices = 0, 
               i = 0)

for (i in seq_along(row.names(mods)))
{
  net2 <- delete_edges(net1, E(net1)[which(E(net1)$weight < q[i])])
  isolated <-  which(degree(net2)==0)
  net2 <-  delete.vertices(net2, isolated)
  mods$mod[i] <- modularity(cluster_infomap(net2))
  mods$edges[i] <- length(E(net2))
  mods$vertices[i] <- length(V(net2))
  mods$i[i] <- i
}

ggplot(mods, aes(i, mod)) +
  geom_point()

ggplot(mods, aes(i, vertices)) +
  geom_point()


net3 <- delete_edges(net1, E(net1)[which(E(net1)$weight < q[45])])
im <- cluster_infomap(net3)
modularity(im)

## EDGE TRIMMING BY "VALUE"

aor <- aor %>% mutate(lor2 = lor^2)
values <- aor %>%  select(lor2) %>% arrange(lor2) %>% filter(lor2 > 0.0) %>% distinct(lor2)

net1 <- graph_from_data_frame(d = aor, vertices = fc, directed = F)
l1 <- layout.fruchterman.reingold(net1)
# V(net1)$size <- 10 * V(net1)$count/max(V(net1)$count)
V(net1)$size <- 1
E(net1)$weight <- aor$net_weights
# E(net1)$width <- 0.5 * abs(aor$lor)
# E(net1)$color <- ifelse(aor$lor > 0, "green", "grey")
plot(net1, vertex.label=NA, layout = l1)


mods <- tibble(mod = rep(0.0, length(row.names(values))),
               edges = 0, 
               vertices = 0, 
               i = 0)

for (i in seq_along(row.names(values)))
{
  net2 <- delete_edges(net1, which(E(net1)$lor2 < values$lor2[i]))
  isolated <-  which(degree(net2)==0)
  net2 <-  delete.vertices(net2, isolated)
  mods$mod[i] <- modularity(cluster_infomap(net2))
  mods$edges[i] <- length(E(net2))
  mods$vertices[i] <- length(V(net2))
  mods$i[i] <- i
}

ggplot(mods, aes(i, mod)) +
  geom_point()

net_low <- delete_edges(net1, E(net1)[which(E(net1)$lor2 < values$lor2[2500])])
im_low <- cluster_infomap(net_low)
modularity(im_low)
plot(im_low, net_low, vertex.label=NA)

net_high <- delete_edges(net1, E(net1)[which(E(net1)$lor2 < values$lor2[2850])])
im_high <- cluster_infomap(net_high)
modularity(im_high)
plot(im_high, net_high, vertex.label=NA)

sink("communities im_low  2020-03-13.txt")
print(communities(im_low))
sink()

net_super <- delete_edges(net1, E(net1)[which(E(net1)$lor2 < values$lor2[3020])])
im_super <- cluster_infomap(net_super)
modularity(im_super)
plot(im_super, net_super, vertex.label=NA)

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


