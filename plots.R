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
                  m2 = paste(".", data_raw[,9], sep=""),
                  l1ctcf = as.logical(data_raw[,10]),
                  m1ctcf = as.logical(data_raw[,11]),
                  l2ctcf = as.logical(data_raw[,12]),
                  m2ctcf = as.logical(data_raw[,13]))
  return(d)
}

# read a file and convert to a dataframe
results2dataframe <- function(file){
  data_raw <- read.table(file, header=F, sep = "\t")
  
  d <- data.frame(QP = data_raw[,1], 
                  l1 = data_raw[,6],
                  m1 = data_raw[,7],
                  l2 = data_raw[,8],
                  m2 = data_raw[,9],
                  l1ctcf = as.logical(data_raw[,10]),
                  m1ctcf = as.logical(data_raw[,11]),
                  l2ctcf = as.logical(data_raw[,12]),
                  m2ctcf = as.logical(data_raw[,13]))
  return(d)
}

# generate k-mer list containing CTCF binding motif
get_ctcf_list <- function(d){
  return(union(union(d$l1[which(d$l1ctcf == TRUE)], 
                     d$l2[which(d$l2ctcf == TRUE)]),
               union(d$m1[which(d$m1ctcf == TRUE)],
                     d$m2[which(d$m2ctcf == TRUE)])))
}

# plot QP weights. horizontal line represents standard deviation
plot_pdf_QP_weights <- function(d, pdfName, out = "None"){
  if(out == "pdf"){
    pdf(pdfName, paper="a4")
  }
  plot(d$QP, xlim = c(0,length(d$QP)), ylim = 1.2 * range(d$QP))
  par(T)
  abline(h = (mean(d$QP) - 2 * sd(d$QP)))
  par(T)
  abline(h = (mean(d$QP) - sd(d$QP)))
  par(T)
  abline(h = (mean(d$QP) + sd(d$QP)))
  par(T)
  abline(h = (mean(d$QP) + 2 * sd(d$QP)))
  par(F)
  if(out == "pdf"){
    dev.off()
  }
}

# select significant k-mer pairs
select_significant <- function(d, x){
  plus  <- which(d$QP <= (mean(d$QP) - x * sd(d$QP)))
  minus <- which(d$QP >= (mean(d$QP) + x * sd(d$QP)))
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
plot_pdf_bi_color <- function(g, ctcf_bi, lay, pdfName, out = "None"){
  # edge color preparation
  cols = colorRampPalette(c("#0068b7","gray","#f39800"))
  # vertex color, shape preparation
  vcols <- c("steelblue", "orange")
  vshapes <- c("circle", "square")
  # write to pdf
  if(out == "pdf"){
    pdf(pdfName, paper="a4")
  }
  # plot with color
  plot(g, layout = lay,
       vertex.size = (4 * (as.numeric(V(g)$label %in% ctcf_bi)) + 1),
       vertex.color = vcols[as.numeric(V(g)$label %in% ctcf_bi) + 1],
       vertex.shape = vshapes[as.numeric(V(g)$label %in% ctcf_bi) + 1],
       edge.color = cols(10)[round(E(g)$weights * 5 + 5)])
  if(out == "pdf"){
    dev.off()
  }
}

# plot a bipartite graph (gray)
plot_pdf_bi_mono <- function(g, ctcf_bi, lay, pdfName, out = "None"){
  # vertex shape preparation
  vshapes <- c("circle", "square")
  # write to pdf
  if(out == "pdf"){
    pdf(pdfName, paper="a4")
  }
  # plot without color
  plot(g, layout = lay,
       vertex.size = (4 * (as.numeric(V(g)$label %in% ctcf_bi)) + 1),
       vertex.color = "gray",
       vertex.shape = vshapes[as.numeric(V(g)$label %in% ctcf_bi) + 1],
       edge.width = abs(E(g)$weights) ^ 10)
  if(out == "pdf"){
    dev.off()
  }
}

# plot a graph (color)
plot_pdf_color <- function(g, ctcf, pdfName, out = "None"){
  # edge color preparation
  cols = colorRampPalette(c("#0068b7","gray","#f39800"))
  # vertex color, shape preparation
  vcols <- c("steelblue", "orange")
  vshapes <- c("circle", "square")
  
  # write to pdf
  if(out == "pdf"){
    pdf(pdfName, paper="a4")
  }
  # plot with color
  plot(g,
       vertex.size = (4 * (as.numeric(V(g)$label %in% ctcf)) + 1),
       vertex.color = vcols[as.numeric(V(g)$label %in% ctcf) + 1],
       vertex.shape = vshapes[as.numeric(V(g)$label %in% ctcf) + 1],
       edge.color = cols(10)[round(E(g)$weights * 5 + 5)])
  if(out == "pdf"){
    dev.off()
  }
}

# plot a graph (gray)
plot_pdf_mono <- function(g, ctcf, pdfName, out = "None"){
  # vertex color, shape preparation
  vshapes <- c("circle", "square")
  
  # write to pdf
  if(out == "pdf"){
    pdf(pdfName, paper="a4")
  }
  # plot without color
  plot(g,
       vertex.size = (4 * (as.numeric(V(g)$label %in% ctcf)) + 1),
       vertex.color = "gray",
       vertex.shape = vshapes[as.numeric(V(g)$label %in% ctcf) + 1],
       edge.width = abs(E(g)$weights) ^ 10)
  if(out == "pdf"){
    dev.off()
  }
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
g_sig_bi <- dataframe2bipartite_graph(select_significant(d_bi, 2))
lay_sig_bi <- computeLayout_bipartite(g_sig_bi)
d <- results2dataframe(file_name)
g <- dataframe2graph(select_significant(d, 2))
ctcf_bi <- get_ctcf_list(d_bi)
ctcf <- get_ctcf_list(d)


############################
# plot to pdf files
############################

out_format = "pdf"

dim(d_bi)
mean(d_bi$QP)
sd(d_bi$QP)
plot_pdf_QP_weights(d_bi, pdf_name_QPweights, out_format)
plot_pdf_bi_color(g_all_bi, ctcf_bi, lay_all_bi, pdf_name_all_bi_col, out_format)
plot_pdf_bi_mono(g_all_bi,  ctcf_bi, lay_all_bi, pdf_name_all_bi_mono, out_format)
plot_pdf_bi_color(g_sig_bi, ctcf_bi, lay_sig_bi, pdf_name_bi_col, out_format)
plot_pdf_bi_mono(g_sig_bi,  ctcf_bi, lay_sig_bi, pdf_name_bi_mono, out_format)
plot_pdf_color(g, ctcf, pdf_name_col, out_format)
plot_pdf_mono(g,  ctcf, pdf_name_mono, out_format)

