library(igraph)

args <- commandArgs(trailingOnly = T)
if(length(args) != 1) {
  write("Error: a commandline argument for <infile> is required", stderr())
  q()
}

############################
# definitions of functions
############################

# read a file and convert to a dataframe
results2dataframe_bipartite <- function(file){
  data_raw <- read.table(file, header=F, sep = "\t")
  
  d <- data.frame(QP = data_raw[,1], 
                  l1 = paste(data_raw[,6], ".", sep=""),
                  m1 = paste(".", data_raw[,7], sep=""),
                  l2 = paste(data_raw[,8], ".", sep=""),
                  m2 = paste(".", data_raw[,9], sep=""))
  return(d)
}

# read a file and convert to a dataframe
results2dataframe <- function(file){
  data_raw <- read.table(file, header=F, sep = "\t")
  
  d <- data.frame(QP = data_raw[,1], 
                  l1 = data_raw[,6],
                  m1 = data_raw[,7],
                  l2 = data_raw[,8],
                  m2 = data_raw[,9])
  return(d)
}

# plot QP weights. horizontal line represents standard deviation
plot_pdf_QP_weights <- function(d, pdfName){
  pdf(pdfName, paper="a4")
  plot(d$QP, xlim = c(0,length(d$QP)), ylim = 1.2 * range(d$QP))
  par(T)
  abline(h = (mean(d$QP) - sd(d$QP)))
  par(T)
  abline(h = (mean(d$QP) + sd(d$QP)))
  par(F)
  dev.off()
}

# select significant k-mer pairs
select_significant <- function(d){
  plus  <- which(d$QP <= (mean(d$QP) - sd(d$QP)))
  minus <- which(d$QP >= (mean(d$QP) + sd(d$QP)))
  return(d[c(plus, minus),])
}

# convert a dataframe into a bipartite graph
dataframe2bipartite_graph <- function(d){
  d1<-data.frame(l = d$l1, m = d$m1, weights = d$QP)
  d2<-data.frame(l = d$l2, m = d$m2, weights = d$QP)
  d_merge <- merge(d1, d2, all = T)
  # convert to a graph
  g <- graph_from_data_frame(d_merge, directed = FALSE)
  
  # normalize weights
  E(g)$weights <- E(g)$weights / max(abs(E(g)$weights)) 
  V(g)$label <- V(g)$name
  
  # bipartite graph
  V(g)[name %in% d_merge$l]$type <- TRUE
  V(g)[name %in% d_merge$m]$type <- FALSE
  
  return(g)
}

# convert a dataframe into a graph
dataframe2graph <- function(d){
  d1<-data.frame(l = d$l1, m = d$m1, weights = d$QP)
  d2<-data.frame(l = d$l2, m = d$m2, weights = d$QP)
  d_merge <- merge(d1, d2, all = T)
  # convert to a graph
  g <- graph_from_data_frame(d_merge, directed = FALSE)
  
  # normalize weights
  E(g)$weights <- E(g)$weights / max(abs(E(g)$weights)) 
  V(g)$label <- V(g)$name
  
  return(g)
}

# compute layout
computeLayout_bipartite <- function(g){
  # graph layout
  lay <- layout_as_bipartite(g)
  lay2 <- lay
  lay2[,1] <- lay[,2]
  lay2[,2] <- lay[,1]
  return(lay2)
}

# plot a bipartite graph (color)
plot_pdf_bi_color <- function(g, lay, pdfName){
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

# plot a bipartite graph (gray)
plot_pdf_bi_mono <- function(g, lay, pdfName){
  # write to pdf
  pdf(pdfName, paper="a4")
  # plot without color
  plot(g, layout = lay,
       vertex.size = 1,
       edge.width = abs(E(g)$weights) ^ 10)
  dev.off()   
}

# plot a graph (color)
plot_pdf_color <- function(g, pdfName){
  # edge color preparation
  cols = colorRampPalette(c("#0068b7","gray","#f39800"))
  # write to pdf
  pdf(pdfName, paper="a4")
  # plot with color
  plot(g,
       vertex.size = 1,
       edge.color = cols(10)[round(E(g)$weights * 5 + 5)])
  dev.off()   
}

# plot a graph (gray)
plot_pdf_mono <- function(g, pdfName){
  # write to pdf
  pdf(pdfName, paper="a4")
  # plot without color
  plot(g,
       vertex.size = 1,
       edge.width = abs(E(g)$weights) ^ 10)
  dev.off()   
}

############################
# file name preparation
############################

file_name <- args[1]
pdf_name_QPweights <- paste(file_name, "dist", "pdf", sep=".")
pdf_name_all_bi_col <- paste(file_name, "all", "bi", "col", "pdf", sep=".")
pdf_name_all_bi_mono <- paste(file_name, "all", "bi", "mono", "pdf", sep=".")
pdf_name_bi_col <- paste(file_name, "bi", "col", "pdf", sep=".")
pdf_name_bi_mono <- paste(file_name, "bi", "mono", "pdf", sep=".")
pdf_name_col <- paste(file_name, "col", "pdf", sep=".")
pdf_name_mono <- paste(file_name, "mono", "pdf", sep=".")

############################
# data preparation
############################

d_bi <- results2dataframe_bipartite(file_name)
g_all_bi <- dataframe2bipartite_graph(d_bi)
lay_all_bi <- computeLayout_bipartite(g_all_bi)
g_sig_bi <- dataframe2bipartite_graph(select_significant(d_bi))
lay_sig_bi <- computeLayout_bipartite(g_sig_bi)
d <- results2dataframe(file_name)
g <- dataframe2graph(select_significant(d))


############################
# plot to pdf files
############################

mean(d_bi$QP)
sd(d_bi$QP)
plot_pdf_QP_weights(d_bi, pdf_name_QPweights)
plot_pdf_bi_color(g_all_bi, lay_all_bi, pdf_name_all_bi_col)
plot_pdf_bi_mono(g_all_bi,  lay_all_bi, pdf_name_all_bi_mono)
plot_pdf_bi_color(g_sig_bi, lay_sig_bi, pdf_name_bi_col)
plot_pdf_bi_mono(g_sig_bi,  lay_sig_bi, pdf_name_bi_mono)
plot_pdf_color(g, pdf_name_col)
plot_pdf_mono(g, pdf_name_mono)

