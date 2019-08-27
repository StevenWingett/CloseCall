

files <- list.files(path=".", pattern="*.cis_distance_by_group.txt", full.names=T, recursive=FALSE)

files

lapply(files, function(file) {
  
  outputfilename=paste(file, "pdf", sep = ".")
  pdf(outputfilename, height=20, width=20) 





  data <- read.delim(file, header=T)
  data <- subset(data, Cis_Distances > 100)
  data <- subset(data, Cis_Distances > 500)
  data <- subset(data, Cis_Distances > 1000)
  

  title <- 'Cis read separation by bacrode group size'
  title <- paste(title, file, sep="\n")
  
  bp <- boxplot(data$Cis_Distances ~ data$Group_Size, xlab="Barcode Group Size", ylab="Separation Distancs (bps)", main=title)
  dev.off()

}
)