#
files <- list.files(path=".", pattern="*.includes_too_large.txt$", full.names=T, recursive=FALSE)

files
lapply(files, function(file) {
  
  title <- 'Frequency of barcodes with n unique reads'
  title <- paste(title, file, sep="\n")
  
  outputfilename=paste(file, "pdf", sep = ".")
  pdf(outputfilename, height=20, width=20)
  
  results <- read.delim(file, header=FALSE, stringsAsFactors=FALSE)
  
  lastbar <- results[1,2]    #Number of allowed bars
  lastbar <- as.numeric(lastbar)
  results <- results[-1,]    #Remove row 1
  
  freq <- table(results[,1])
  freq <- as.matrix(freq)
  freq <- table(freq[,1])
  
  max.value = max(freq)
  y.graph = 1.2 * max.value
  
  bp <-barplot(freq, main=title, xlab='n', ylab='Frequency', ylim=c(0, y.graph) )
  
  #Determine position for ab line (if needed)
  if(lastbar < nrow(bp)){
     cutoffline <- (bp[lastbar,1] + bp[(lastbar + 1), 1])/2 
     abline(v=cutoffline)
  }
     text (bp, freq, freq, cex=1, pos=3)
     par(oma=c(0,0,0,0))
     
     dev.off()
}

 )