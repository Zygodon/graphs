
# Libraries
library("RMySQL")
library(tidyverse)
library(networkD3)
library(igraph)
library(plotly)
library(ggpmisc)
library(DT)
library(gridExtra)

pwor <- read.csv("pwor_10K.csv", header = TRUE, sep = ",")


friends <- (pwor %>% filter((sim_q_high < ci_low.x) & (sfx == "yes"))
            %>% select(A, B, qor, share_2x2)
            %>% rename("log(odds ratio)" = qor) 
            %>% rename("shared quadrats" = share_2x2))



fp <- friends %>% select(A, B) %>% pivot_longer(cols = c(A,B), names_to = "origin", values_to = "species")
fc <- fp %>% group_by(species) %>% summarise(count = n())

net1 <- graph_from_data_frame(d = friends, vertices = fc, directed = F)
l <- layout.fruchterman.reingold(net1)
plot(net1, layout=l, vertex.size = 2, vertex.label=NA)
V(net1)$size <- V(net1)$count
plot(net1, vertex.label=NA)


wc <- cluster_walktrap(net1)
modularity(wc)
plot(wc, net1, vertex.label=NA)
