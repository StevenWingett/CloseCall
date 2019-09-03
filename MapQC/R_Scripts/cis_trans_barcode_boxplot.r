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