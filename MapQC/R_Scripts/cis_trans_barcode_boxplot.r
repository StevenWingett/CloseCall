files <- list.files(path=".", pattern="*.cis_trans_by_group.txt", full.names=T, recursive=FALSE)

lapply(files, function(file) {
  
  
  
  
  graph.title = paste("%Trans Vs Barcode Group Size", file, sep="\n")
  outputfilename=paste(file, "pdf", sep = ".")
  pdf(outputfilename, height=20, width=20) 
  
  
  data <- read.delim(file)
  boxplot(data$Percentage_Trans ~ data$Group_Size, main=graph.title,
          xlab="Barcode Group Size", ylab="%Trans"
  )
  
  
  
  
  dev.off()
  
}
)