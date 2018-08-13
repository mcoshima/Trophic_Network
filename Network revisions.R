setwd("C:/Users/w986430/OneDrive - The University of Southern Mississippi/WD/Manuscript Figs/Manuscript Revisions")
library(igraph)
library(NetIndices)
library(keyplayer)
library(tibble)
library(enaR)

prelim.network.2 <- f.sp

#Predator species name, lowest.taxa, Freq.of.Occur., Dry Weight, Weight, Volume, IRI lowest taxa level
el <- matrix(NA, nrow = dim(prelim.network.2)[1], ncol = 3)
el[,1] <- as.character(prelim.network.2[,1])
el[,2] <- as.character(prelim.network.2[,2])
el[,3] <-  as.character(prelim.network.2[,3])


g.2 <- graph.edgelist(el[,1:2])
E(g.2)$weight <- as.numeric(el[,3]) 
g.2 <- simplify(g.2, remove.multiple = T )

 # org.vals <- matrix(data = NA, nrow = 8, ncol = 9)
 # colnames(org.vals) <- c("robust", "con.w", "LD.w", "cc", "cc.w", "apl.w", "con", "ld", "S")
 # rownames(org.vals) <- c("Species.v", "Genus.v", "Family.v", "Order.v", "Order.w", "Family.w","Genus.w","Species.w")
 # g.2.adj.mat <- as_adjacency_matrix(g.2, attr = "weight")
 # #robustness, weighted connectance and weighted link density
 # org.vals[a,c(1:3)] <- oascend.funs(g.2.adj.mat)
 # #clustering coefficient tells how much a sub-group of nodes are connected
 # org.vals[a,4] <- transitivity(g.2, type = "average")
 # #weighted cc
 # g.2.undir <- as.undirected(g.2, mode = "collapse", edge.attr.comb = "sum")
 # org.vals[a,5] <- mean(transitivity(g.2.undir, type = "weighted"), na.rm = TRUE)
 # #weighted average path length
 # org.vals[a,6] <- oapl.fun(g.2)
 # #connectance - the proportion of possible links between species that are realized (L/S^2)
 # org.vals[a,7] <- edge_density(g.2)
 # #link density
 # org.vals[a,8] <- length(E(g.2))/length(V(g.2))
 # org.vals[a,9] <- length(V(g.2))

# deg <- rownames_to_column(as.data.frame(degree(g.2, v=V(g.2))))
# deg.ranks <- matrix(data = NA, nrow = nrow(deg), ncol = 10)
# for (j in 1:10){
# deg$rank <- rank(deg$`degree(g.2, v = V(g.2))`, ties.method = "random")
# deg <- deg[order(-deg$rank),]
# deg.ranks[,j] <- deg$rowname
# }

# close <- rownames_to_column(as.data.frame(closeness(g.2, vids = V(g.2), mode = "total", normalized = TRUE)))
# close.ranks <- matrix(data = NA, nrow = nrow(close), ncol = 10)
# for (j in 1:10){
# close$rank <- rank(round(close$`closeness(g.2, vids = V(g.2), mode = "total", normalized = TRUE)`,2), ties.method = "random")
# close <- close[order(close$rank),]
# close.ranks[,j] <- close$rowname
# }

# vert <- as.vector(V(g.2)$name)
# Frag.df <- data.frame("Node" = vert, "Frag" = fragment(adj.matrix = g.2.adj.mat))
# frag.ranks <- matrix(data = NA, nrow = nrow(Frag.df), ncol = 10)
# for (j in 1:10){
#   Frag.df$rank <- Frag.df$rank <- rank(round(Frag.df$fragment,3), ties.method = "random")
#   Frag.df <- Frag.df[order(-Frag.df$rank),]
#   Frag.df$Node <- as.character(Frag.df$Node)
#   frag.ranks[,j] <- Frag.df$Node
# }

# wd <- rownames_to_column(as.data.frame(strength(g.2, vids = V(g.2))))
# wd.ranks <- matrix(data = NA, nrow = nrow(wd), ncol = 10)
# for (j in 1:10){
#   wd$rank <- rank(round(wd$`strength(g.2, vids = V(g.2))`,digits = 0), ties.method = "random")
#   wd <- wd[order(-wd$rank),]
#   wd.ranks[,j] <- wd$rowname
#   }


# node.ranks <- matrix(data = NA, nrow = gorder(g.2), ncol = 10)
# v.vec <- V(g.2)$name
# for(j in 1:10){
#   order.node <- sample(v.vec, size = length(v.vec), replace = F)
#   node.ranks[,j] <- order.node
# }

node.removed.vec <- c()
ind.vals <- matrix(data = NA, nrow = gorder(g.2), ncol = 9)
colnames(ind.vals) <- c("robust", "con.w", "LD.w", "cc", "cc.w", "apl.w", "con", "ld", "f")
removals <- list()

g.1 <- g.2
for (i in 1:nrow(deg.ranks)){
  node.remove <- V(g.1)[deg.ranks[i,j]]#find node with most connections (A)
  node.removed.vec[i] <- node.remove
     #get edges to delete, edges that are pointing in towards the node
    if(length(unlist((incident_edges(g.1, node.remove, mode = "in"))))==0){
      g.9 = delete.vertices(g.1, node.remove)
      g.9.adj.mat <- as_adjacency_matrix(g.9, attr = "weight")
      #robustness, weighted connectance and weighted link density
      ind.vals[i,] <- ascend.funs(g.9.adj.mat)
      #clustering coefficient tells how much a sub-group of nodes are connected
      ind.vals[i,4] <- transitivity(g.9, type = "average")
      #weighted cc
      g.9.undir <- as.undirected(g.9, mode = "collapse", edge.attr.comb = "sum")
      ind.vals[i,5] <- mean(transitivity(g.9.undir, type = "weighted"), na.rm = TRUE)
      #weighted average path length
      ind.vals[i,6] <- apl.fun(g.9)
      #connectance - the proportion of possible links between species that are realized (L/S^2)
      ind.vals[i,7] <- edge_density(g.9)
      #link density
      ind.vals[i,8] <- length(E(g.9))/length(V(g.9))
      #fraction of nodes removed
      ind.vals[i,9] <- (length(V(g.2))-length(V(g.9)))/length(V(g.2))
      g.1 <- g.9
    } else{
    edges.in.delete <-
      unlist(incident_edges(g.1, node.remove, mode = "in"))
    #get weights of edges that will be deleted
    edge.weights <-
      edge_attr(g.1, "weight", index = E(g.1)[edges.in.delete])
    
    for (m in 1:length(edges.in.delete)) {
      #selects the edges that will be deleted after removing node
      delete.this.edge <- edges.in.delete[m]
      #selects neighbors that are feeding on A
      removed.node.neigh <-
        unlist(ego(g.1, 1, nodes = node.remove, mode = "in"))
      #removed.node.neigh[1] is node that is being removed (A)
      removed.node.neigh <- removed.node.neigh[!(removed.node.neigh == node.remove)]
      
      #selects the edges of the removed node's neighbors
      cousin.edge.ind <-
        unlist(incident_edges(g.1, removed.node.neigh[m], mode = "out"))
      cousin.edge.ind <-
        cousin.edge.ind[!(cousin.edge.ind == delete.this.edge)]
      #sum of those edge weights
      edge.sum <- sum(E(g.1)[cousin.edge.ind]$weight)
      #determines the weight to redistribute to each edge from ratio
      weight.to.add <-
        E(g.1)[delete.this.edge]$weight * (E(g.1)[cousin.edge.ind]$weight / edge.sum) 
      #calculates new weight
      new.weight <- E(g.1)[cousin.edge.ind]$weight + weight.to.add 
      #assigns new weights to edges
      E(g.1)[cousin.edge.ind]$weight <- new.weight
      
    }
  
  
  #removes node
    g.9 = delete.vertices(g.1, node.remove)
    g.9.adj.mat <- as_adjacency_matrix(g.9, attr = "weight")
    #robustness, weighted connectance and weighted link density
    ind.vals[i,] <- ascend.funs(g.9.adj.mat)
    #clustering coefficient tells how much a sub-group of nodes are connected
    ind.vals[i,4] <- transitivity(g.9, type = "average")
    #weighted cc
    g.9.undir <- as.undirected(g.9, mode = "collapse", edge.attr.comb = "sum")
    ind.vals[i,5] <- mean(transitivity(g.9.undir, type = "weighted"), na.rm = TRUE)
    #weighted average path length
    ind.vals[i,6] <- apl.fun(g.9)
    #connectance - the proportion of possible links between species that are realized (L/S^2)
    ind.vals[i,7] <- edge_density(g.9)
    #link density
    ind.vals[i,8] <- length(E(g.9))/length(V(g.9))
    #fraction of nodes removed
    ind.vals[i,9] <- (length(V(g.2))-length(V(g.9)))/length(V(g.2))
    g.1 <- g.9
}
}


removals[[j]] <- ind.vals

ave.mat <- list()
for(m in 1:ncol(removals[[1]])){
  
  ave.mat[[m]] <- list.averages.fun(removals, m)
}
means <- do.call(cbind, ave.mat)
colnames(means) <- c("robust.m", "robust.sd","con.w.m", "con.w.sd","LD.w.m", "LD.w.sd","cc.m", "cc.sd","cc.w.m", "cc.w.sd","apl.w.m", "apl.w.sd","con.m","con.sd", "ld.m","ld.sd", "f.m", "f.sd")


write.csv(means, file = "weight.deg.means.or.csv", row.names = F)
write.csv(means, file = "fo.close.means.or.csv", row.names = F)
write.csv(means, file = "fo.frag.means.or.csv", row.names = F)
write.csv(means, file = "fo.wd.means.or.csv", row.names = F)

write.csv(org.vals, file = "vw.original.vals.csv")
