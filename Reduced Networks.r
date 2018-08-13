install.packages("ritis")
install.packages("igraph")
library(igraph)
library(dplyr)
library(ritis)
## Function to wrap long strings
# Source: http://stackoverflow.com/a/7367534/496488
wrap_strings <- function(vector_of_strings,width){
  as.character(sapply(vector_of_strings, FUN=function(x){
    paste(strwrap(x, width=width), collapse="\n")
  }))
}

vl <- read.csv(file.choose())
# Apply the function to wrap the node labels
vl = wrap_strings(vl, 8)


# Function to increase node separation (for explanatory details, see the link below)
# Source: http://stackoverflow.com/a/28722680/496488
layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
  E(g)$weight <- 1
  
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
  
  l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
  return(l)
}



####For making manuscript network figure
setwd("C:/Users/w986430/OneDrive - The University of Southern Mississippi/WD/Manuscript Figs/Manuscript Revisions/Network Revisions R")
fo.sp <- read.csv("fo.sp.csv", stringsAsFactors = F)
head(fo.sp)
#fo.tax <- read.csv(".dfo.tax.fixed.csv", stringsAsFactors = F)
fo.tax <- read.csv(file.choose(), stringsAsFactors = F)

for(i in 1:nrow(fo.sp)){
  pred.row <- which(fo.sp$Species[i] == fo.tax$Species)
  pred.row <- pred.row[1]
  fo.sp[i,1] <- paste(strsplit(fo.tax[pred.row,5], split = "")[[1]][1], ". ", fo.sp[i,1],  sep = "")
}

for(i in 1:nrow(fo.sp)){
  prey.row <- which(fo.sp$Species.2[i] == fo.tax$Species.2)
  prey.row <- prey.row[1]
  fo.sp[i,2] <- paste(strsplit(fo.tax$Genus.1[prey.row], split = "")[[1]][1], ". ", fo.sp[i,2],  sep = "")
}

el.mat.2 <- as.matrix(fo.sp)
g.2 <- graph.edgelist(el.mat.2[,1:2])
E(g.2)$weight <- as.numeric(el.mat.2[,3]) #adds weight
g.2 <- simplify(g.2, remove.multiple = T )

men.pred <- adjacent_vertices(g.2, v = V(g.2)[24], mode = "all")
men.pred <- unlist(men.pred)

pat.pred <- get.edges(g.2, es = E(g.2)[from(men.pred)])
pat.subg2 <- ends(g.2, es = E(g.2)[from(pat.pred)])
vps <- array(t(pat.subg2))
pat.subg2 <- cbind(pat.subg2,  E(g.2)[from(pat.pred)]$weight)
g2.sub <- graph.edgelist(pat.subg2[,1:2])
E(g2.sub)$weight <- as.numeric(pat.subg2[,3]) #adds weight
g2.sub <- simplify(g2.sub, remove.multiple = T )
e.list <- get.edge.ids(g.2, vp = vps)
e.weights <- E(g2.sub)$weight

vsize <- log(strength(g2.sub, vids = V(g2.sub)))
l.1 <- layout.by.attr(g2.sub, wc = 15)
move.node.labels <- V(g2.sub)$name[c(52,49)]

label.dist <- rep(0,65)
label.dist[c(10,52,49)] <- .5

# communities <- cluster_optimal(g2.sub)
# membership(communities)
# 
# col. <- c("#FFC9B5", "#493548", "#99D19C", "#79C7C5", "#B2B1CF", "#474448", "#ED9D61")
# 
# col. <- rainbow(10)
# tiff(file="manuscript.net.fig.tiff" , height=10,  width=8, pointsize=12 , units = "in", res = 500)
# par(mar = rep(1,4), oma = rep(1.5,4))
# plot(g2.sub, vertex.size=vsize*3,     
#      vertex.color= col.[membership(communities)],
#      vertex.label = ifelse(vsize > 3, V(g2.sub)$name, NA),
#      vertex.shape = "sphere",
#      #vertex.label.color = "black",
#      vertex.label.family = "serif",
#      vertex.label.font = 3,
#      #vertex.label.dist =label.dist,
#      edge.width = log(e.weights)*2,
#      edge.arrow.size = 0.5,
#      edge.color = "gray",
#      edge.curved = T,
#      layout=l.1)
# dev.off()

tiff(file="manuscript.net.fig.tiff" , height=10,  width=8, pointsize=12 , units = "in", res = 500)
par(mar = rep(1,4), oma = rep(1.5,4))
plot(g2.sub, vertex.size=as.numeric(vsize),     
     vertex.color= "white",
     vertex.label = V(g2.sub)$name,
     vertex.label.color = "black",
     vertex.label.family = "serif",
     vertex.label.font = 3,
     vertex.label.dist =label.dist,
     edge.width = log(e.weights)*2,
     edge.arrow.size = 0.5,
     edge.color = "gray",
     layout=l.1)
dev.off()

example.adj.mat <- as_adjacency_matrix(g2.sub, attr = "weight")
example.adj.mat <- as.matrix(example.adj.mat)
example.adj.mat <- as.data.frame(example.adj.mat)
write.csv(example.adj.mat, file = "manuscritp.example.adj.matrix.csv", row.names = T)

##Apendix table
head(fo.tax)
fo.sp <- read.csv("fo.sp.csv", stringsAsFactors = F)
pred.vec <- pred.order <- pred.class <- c()
for(i in 1:nrow(fo.sp)){
  pred.row <- which(fo.sp$Species[i] == fo.tax$Species)
  pred.row <- pred.row[1]
  pred.vec[i] <- fo.tax[pred.row,1]
  pred.order[i] <- fo.tax[pred.row,3]
  pred.class[i] <- fo.tax[pred.row,2]
}

prey.vec <- prey.order <- prey.class <- c()
for(i in 1:nrow(fo.sp)){
  prey.row <- which(fo.sp$Species.2[i] == fo.tax$Species.2)
  prey.row <- prey.row[1]
  prey.vec[i] <- fo.tax[prey.row,12]
  prey.order[i] <- fo.tax[prey.row,8]
  prey.class[i] <- fo.tax[prey.row,7]
}


pred.df <- cbind(pred.class, pred.order, pred.vec)
prey.df <- cbind(prey.class, prey.order, prey.vec)
uni.pred <- distinct(as.data.frame(pred.df))
uni.prey <- distinct(as.data.frame(prey.df))
colnames(uni.pred) <- c("class", "order", "species")
colnames(uni.prey) <- c("class", "order", "species")
uni.species <- bind_rows(uni.pred, uni.prey)
#order for Cylinchella bidnetata
uni.species[117,2] <- "Cephalaspidea"
uni.species <- uni.species %>% group_by(order)


library(rvest)
library(dplyr)
library(stringr)
library(ritis)

for(i in 1:nrow(uni.species)){
  spec. <- gsub(pattern =  " ", "-", uni.species[i,3]) 
  uni.species[i,4] <- paste("http://www.fishbase.org/summary/", spec., ".html", sep = "")
}
uni.species[15,4] <- "http://www.fishbase.org/summary/Micropogonias-undulatus.html"  
uni.species[14,4] <- "http://www.fishbase.org/summary/Malacosteus-niger.html"
uni.species <- as.data.frame(uni.species)
uni.species <- uni.species[-102,]

full.df <- matrix(NA, nrow = nrow(uni.species), ncol = 2)
for (i in 1:nrow(uni.species)){
  #created a web ID with the URL from fish.df$URL
  web.id <- as.character(uni.species[i,"V4"]) 
  #opened the html code for the webpage 
  web.html <- read_html(web.id) 
  #pulled information from the table on the webpage
  species.id <- as.character(html_text(html_nodes(web.html, "#ss-sciname")))
  #split after ID= and called the number (2nd object)
  if(length(species.id) == 0){
    next}
  species.id <- unlist(strsplit(species.id, "\t"))
  species.id <- species.id[-which(species.id == "")]
  species.id <- species.id[-which(species.id == "\r\n")]
  full.df[i,] <- species.id
  print(i)
}


miss.auth <- uni.species[which(is.na(full.df[,1])),2]
miss.row <- which(is.na(full.df[,1]))

for (i in 1:length(miss.auth)){
  tax.df <- as.data.frame(search_scientific(miss.auth[i], wt = "json", raw = F))
  if(nrow(tax.df)==0){
    next
  }
  auth. <- as.character(tax.df[1,1])
  tsn. <- tax.df[1,4]
  full.df[miss.row[i],1] <- auth.
  full.df[miss.row[i],2] <- as.character(common_names(tsn = tsn., wt = "json", raw = F)[1,1])
  print(i)
}

sci.name <- uni.species[,3]
tsn. <- common.name <- c()
for (i in 1:length(sci.name)){
  tax.df <- as.data.frame(search_scientific(sci.name[i], wt = "json", raw = F))
  if(nrow(tax.df)==0){
    next
  }
  tsn.[i] <- tax.df[1,4]
  common.name[i] <- as.character(common_names(tsn = tsn.[i], wt = "json", raw = F))
  print(i)
}







#capitallize first letter in the name
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
full.df[,2] <- firstup(full.df[,2])
append2 <- cbind(uni.species, full.df)
append2[which(append2$species == "Plesionika richardi"),3] <- "Stylopandalus richardi"
write.csv(append2, file = "appendix2.network.manuscript.csv", row.names = F)
####

n.vsize = (strength(g.9, vids = V(g.9)) - strength(g2.sub, vids = V(g2.sub)))/strength(g2.sub, vids = V(g2.sub))
nweight <- log(E(g.9)$weight)

png(file="bc.red.net.png" , height=12,  width=13, pointsize=18 , units = "in", res = 800)
par(mar = c(0,0,0,0), oma = rep(0,4), bg = NA, fg = "white")
plot(g.9, 
     vertex.size = ifelse(n.vsize < 0, 1, n.vsize*20),     
     vertex.color= node.col,
     vertex.label = ifelse(n.vsize > .5, vl, NA),
     vertex.label.color = "#EA8425",
     vertex.label.family = "sans",
     vertex.label.font = 2,
     vertex.label.cex = 1.5,
     edge.width = ifelse(nweight < 1, .5, 2*nweight),
     edge.arrow.size = 0.5,
     edge.color = edge.color,
     layout=l.1)
dev.off()

m.vsize = (strength(g.1, vids = V(g.1)) - strength(g2.sub, vids = V(g2.sub)))/strength(g2.sub, vids = V(g2.sub))
mweight <- log(E(g.1)$weight)


png(file="gm.red.net.png" , height=12,  width=13, pointsize=18 , units = "in", res = 800)
par(mar = c(0,0,0,0), oma = rep(0,4), bg = NA, fg = "white")
plot(g.1, 
     vertex.size = ifelse(m.vsize < 0, 1, m.vsize*10),     
     vertex.color= node.col,
     vertex.label = ifelse(m.vsize > .5, vl, NA),
     vertex.label.color = "#EA8425",
     vertex.label.family = "sans",
     vertex.label.font = 2,
     vertex.label.cex = ifelse(m.vsize > .8, 1.35, 1),
     edge.width = ifelse(mweight < 1, .5, 2*mweight),
     edge.arrow.size = 0.5,
     edge.color = edge.color,
     layout=l.1)
dev.off()

