###################################################################################
###################################################################################
##This file is Copyright (C) 2019, Steven Wingett (steven.wingett@babraham.ac.uk)##
##                                                                               ##
##                                                                               ##
##This file is part of CloseCall.                                                ##
##                                                                               ##
##CloseCall is free software: you can redistribute it and/or modify              ##
##it under the terms of the GNU General Public License as published by           ##
##the Free Software Foundation, either version 3 of the License, or              ##
##(at your option) any later version.                                            ##
##                                                                               ##
##CloseCall is distributed in the hope that it will be useful,                   ##
##but WITHOUT ANY WARRANTY; without even the implied warranty of                 ##
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  ##
##GNU General Public License for more details.                                   ##
##                                                                               ##
##You should have received a copy of the GNU General Public License              ##
##along with CloseCall.  If not, see <http://www.gnu.org/licenses/>.             ##
###################################################################################
###################################################################################


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