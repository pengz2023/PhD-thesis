# Draw gene regulatory networks
library(Matrix)
# install.packages("igraph")
library(igraph)
draw_network <- function(adjMatrix){
  graph <- graph_from_adjacency_matrix(adjMatrix, mode = "undirected", diag = FALSE )
  E(graph)$width <- 1
  index <- sort(degree(graph), index.return=TRUE, decreasing = TRUE)
  # top 5 most highly connected nodes ("hub nodes")
  index.sig <- index$ix[1:5]
  index.nonsig <- index$ix[-(1:5)]
  # nodes without edges
  index.nonedges <- index$ix[index$x==0] 
  # hide the indices of nonsignificant nodes
  # V(graph)$name[index.nonsig] <- NA
  # V(graph)$color <- degree(graph)
  V(graph)$color[index.sig] <- "red"
  V(graph)$color[index.nonsig] <- "orange"
  V(graph)$label.cex = 0.8
  V(graph)$label.color = "black"
  # graph <- delete_vertices(graph, index.nonedges)
  node.size <- 1.5*degree(graph)
  plot(graph, vertex.size = node.size, layout = layout.fruchterman.reingold)
}

# network reconstruction for GNBP
path <- './adjacenyMatrix/'
adjMatrix_GNBP <- as.matrix(read.csv(file=paste0(path, "GNBP_adjacenyMatrix.csv"), header=F))
draw_network(adjMatrix_GNBP)

# network reconstruction for tGNBP
path <- './adjacenyMatrix/'
adjMatrix_tGNBP <- as.matrix(read.csv(file=paste0(path, "tGNBP_adjacenyMatrix.csv"), header=F))
draw_network(adjMatrix_tGNBP)


# network reconstruction for quasiGNBP
path <- './adjacenyMatrix/'
adjMatrix_quasiGNBP <- as.matrix(read.csv(file=paste0(path, "quasiGNBP_adjacenyMatrix.csv"), header=F))
draw_network(adjMatrix_quasiGNBP)


# network reconstruction for quasiGNBP-diag
path <- './adjacenyMatrix/'
adjMatrix_quasiGNBP_diag <- as.matrix(read.csv(file=paste0(path, "quasiGNBP_diag_adjacenyMatrix.csv"), header=F))
draw_network(adjMatrix_quasiGNBP_diag)


# network reconstruction for quasi-tGNBP-diag
path <- './adjacenyMatrix/'
adjMatrix_quasi_tGNBP_diag <- as.matrix(read.csv(file=paste0(path, "quasi_tGNBP_diag_adjacenyMatrix.csv"), header=F))
draw_network(adjMatrix_quasi_tGNBP_diag)









