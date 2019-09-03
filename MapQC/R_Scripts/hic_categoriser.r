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


library("RColorBrewer")


files <- list.files(path=".", pattern="*.hic_category_v_bargroup.txt$", full.names=T, recursive=FALSE)

#file <- 'lane3130_JM_Mar2015_fr.4_cDNA_NoIndex_L001_R1.fastq.gz.present.fastq.gz.barcoded.fastq_trim.txt.sam.data_table_include_non_interacting.txt.hic_category_v_bargroup.txt'
freq = ''

lapply(files, function(file) {
  

  
  outputfilename=paste(file, "pdf", sep = ".")
  # pdf(outputfilename, height=20, width=20)
  pdf(outputfilename)
  title <- '%Virtual Di-tag Categories by Bacrode Group Size'
  title <- paste(title, file, sep="\n")
  
  data <- read.delim(file, header=T)
  freq <- table(data)
  freq <- t(freq)
  
  nos <- as.numeric(colnames(freq)[1:ncol(freq)-1])
  ordered.indices <- c(order(nos),ncol(freq))
  ordered.freq <- freq[,ordered.indices]
  freq <- ordered.freq
  
  freq <- apply(freq, 2, function(x){(x/sum(x))*100})
  
  
  layout( matrix(c(1,2), nrow=1), widths = c(0.7, 0.3)   )
  par( mar = c(5,4,5,0) , oma=c(0,0,0,0)  )
  
  
  
  barplot(freq, legend=F, col=brewer.pal(7,"Set1"), xlab="Barcode Group Size (N)",
          ylab = '%', main=title, cex.names=0.8, cex.main=0.7 )
  
  #par(new=TRUE)
  par(oma=c(0,0,0,0), mar=c(0,0,0,1), new=T)
  #Plot an invisible graph here so the legend is visible and a nice size
  plot(0, 0, type='n', bty='n', xaxt="n", yaxt="n", col="transparent")
  
  
  legend<-rev(rownames(freq))
  legend("topleft", y=NULL, inset=c(-0.01,0.15), legend=legend, bty="n", fill=rev(brewer.pal(7,"Set1") ) )
  #
  
  dev.off()
  
}
)


