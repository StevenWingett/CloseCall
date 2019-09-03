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