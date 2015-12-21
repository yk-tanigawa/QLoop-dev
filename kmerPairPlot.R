library(igraph)

args <- commandArgs(trailingOnly = T)
if(length(args) != 1) {
  write("Error: a commandline argument for <infile> is required", stderr())
  q()
}

#################################################
# definition of several functions
AdaBoost2dataframe <- function(file){
  data_raw <- read.table(file, header=F, sep = "\t")

  d <- data.frame(kmer1 = paste(data_raw[,2], ":", sep=""),
                  kmer2 = paste(":", data_raw[,3], sep=""),
                  weights = (-1) ^ data_raw[,5] * data_raw[,1])
  return(d)
}
dataframe2graph <- function(d){
  # convert to a graph
  g <- graph_from_data_frame(d, directed = FALSE)
  
  # normalize weights
  E(g)$weights <- E(g)$weights / max(abs(E(g)$weights)) 
  
  V(g)$label <- V(g)$name
  
  # bipartite graph
  V(g)[name %in% d$kmer1]$type <- TRUE
  V(g)[name %in% d$kmer2]$type <- FALSE
  
  return(g)
}
computeLayout <- function(g){
  # graph layout
  lay <- layout_as_bipartite(g)
  lay2 <- lay
  lay2[,1] <- lay[,2]
  lay2[,2] <- lay[,1]
  return(lay2)
}
plotColor <- function(g, lay, pdfName){
  # edge color preparation
  cols = colorRampPalette(c("#0068b7","gray","#f39800"))
  # write to pdf
  pdf(pdfName, paper="a4")
  # plot with color
  plot(g, layout = lay,
       vertex.size = 1,
       edge.color = cols(10)[round(E(g)$weights * 5 + 5)])
  dev.off()   
}
plotMono <- function(g, lay, pdfName){
  # write to pdf
  pdf(pdfName, paper="a4")
  # plot without color
  plot(g, layout = lay,
       vertex.size = 1,
       edge.width = abs(E(g)$weights) ^ 10)
  dev.off()   
}
#################################################


#fileName <- "../data/chr21.m10000.M100000.KR.KR.dat.k4.t3.886288.T100.stamps"
fileName <- args[1]

d <- AdaBoost2dataframe(fileName)
g <- dataframe2graph(d)
lay <- computeLayout(g)

plotColor(g, lay, paste(fileName, "col", "pdf", sep = "."))
plotMono(g, lay, paste(fileName, "mono", "pdf", sep = "."))
