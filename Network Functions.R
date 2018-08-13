oapl.fun<- function(gobj.){
  obj. <- gobj.
  apl.mat<-matrix(data = NA, nrow=gorder(obj.), ncol = gorder(obj.))
  E(obj.)$weight<- E(obj.)$weight/mean(E(obj.)$weight) #normalize weights, allows comparison across networks with different ranges of tie weights
  for (y in 1:gorder(obj.)){
    apl.mat[y,] <- distances(obj., v = V(obj.)[y], to = V(obj.),  mode = "all", weights = NULL, algorithm = "dijkstra")
    apl.<-mean(apl.mat)
  }
  print(apl.)
}

apl.fun<- function(gobj.){
  obj. <- gobj.
  apl.mat<-matrix(data = NA, nrow=gorder(obj.), ncol = gorder(obj.))
  E(obj.)$weight<- E(obj.)$weight/mean(E(obj.)$weight) #normalize weights, allows comparison across networks with different ranges of tie weights
  for (y in 1:gorder(obj.)){
    apl.mat[y,] <- distances(obj., v = V(obj.)[y], to = V(obj.),  mode = "all", weights = NULL, algorithm = "dijkstra")
    apl.<-mean(apl.mat[-which(is.infinite(apl.mat))])
  }
  print(apl.)
}

ascend.funs <- function(g.9.adj.mat){
  net.obj.9 <- suppressWarnings(pack(flow = g.9.adj.mat))
  ascend <- enaAscendency(net.obj.9)
  robust <- ascend[1,9]
  con.w <- exp((ascend[1,1] - ascend[1,2])/2)/gorder(g.9)
  ld.w <- exp((ascend[1,1] - ascend[1,2])/2)
  ind. <- c(robust, con.w, ld.w)
  print(ind.)
}


oascend.funs <- function(gadj.mat){
  net.obj.9 <- suppressWarnings(pack(flow = gadj.mat))
  ascend <- enaAscendency(net.obj.9)
  robust <- ascend[1,9]
  con.w <- exp((ascend[1,1] - ascend[1,2])/2)/gorder(g.9)
  ld.w <- exp((ascend[1,1] - ascend[1,2])/2)
  ind. <- c(robust, con.w, ld.w)
  print(ind.)
}


list.averages.fun <- function(list.,n){
  list.ave.mat <- matrix(data = NA, nrow = nrow(list.[[1]]), ncol = length(list.))
  for(x in 1:length(list.)){
    list.ave.mat[,x] <- list.[[x]][,n]
    }
  ave <- rowMeans(list.ave.mat, na.rm = T)
  sd. <- apply(list.ave.mat, 1, sd, na.rm = T)
  cbind(ave, sd.)
}


