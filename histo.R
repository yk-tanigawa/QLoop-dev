args <- commandArgs(trailingOnly = T)
if(length(args) != 1) {
  write("Error: a commandline argument for <infile> is required", stderr())
  q()
}

histo_plot_pdf <- function(file_name, pdf_name){
  histo <- read.table(file_name, header=F)
  names(histo) <- c("Permil", "Interaction Intensities")
  pdf(pdf_name, paper="a4")
  plot(histo, log="y")
  dev.off()
}

file_name <- args[1]
pdf_name <- paste(file_name, "pdf", sep=".")
histo_plot_pdf(file_name, pdf_name)


