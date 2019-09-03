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


###################
#Features count all
files <- list.files(path=".", pattern="\\.MonteCarloFeaturesCount\\.", full.names=T, recursive=FALSE)

lapply(files, function(file) {
  display <- paste("Processing", file, sep=" ")
  print(display)
  results <- read.delim(file, header=TRUE, stringsAsFactors=FALSE)
  
  outputfilename=paste(file, "pdf", sep = ".")
  pdf(outputfilename, height=20, width=20)
  
  xmax <- log( (max(results$Simulated_Average_Frequency) * 1.2), 10)
  ymax <- log( max(results$Observed_Frequency * 1.2), 10)
  title <- paste("Feature Count Scatter Plot (observed vs simulated):", file, sep=" ")
  plot( cbind(  log(results$Simulated_Average_Frequency, 10), log(results$Observed_Frequency, 10) ),
        xlim=c(0, xmax), ylim=c(0, ymax),   xaxs="i",yaxs="i", xlab="Log10(simulated)", ylab="Log10(observed)", 
        main=title, cex.lab=1.5, cex.axis=1.5)
  abline(0,1)
  cex.main=3
  
 dev.off()
})



###################
#Features count Valency greater than 1
files <- list.files(path=".", pattern="\\.MonteCarloFeaturesCountValGt1\\.", full.names=T, recursive=FALSE)

lapply(files, function(file) {
  display <- paste("Processing", file, sep=" ")
  print(display)
  results <- read.delim(file, header=TRUE, stringsAsFactors=FALSE)
  
  outputfilename=paste(file, "pdf", sep = ".")
  pdf(outputfilename, height=20, width=20)
  
  xmax <- log( (max(results$Simulated_Average_Frequency) * 1.2), 10)
  ymax <- log( max(results$Observed_Frequency * 1.2), 10)
  title <- paste("Interacting Features Count Scatter Plot (observed vs simulated):", file, sep="\n")
  plot( cbind(  log(results$Simulated_Average_Frequency, 10), log(results$Observed_Frequency, 10) ),
        xlim=c(0, xmax), ylim=c(0, ymax),   xaxs="i",yaxs="i", xlab="Log10(simulated)", ylab="Log10(observed)", 
        main=title, cex.lab=1.5, cex.axis=1.5)
  abline(0,1)
  cex.main=3
  
  dev.off()
})



#####################################
#Unallocated pool results
files <- list.files(path=".", pattern="\\.MonteCarloUnallocatedPools\\.", full.names=T, recursive=FALSE)
loop <<- 1
graph.resuls <<- NULL

lapply(files, function(file) {
  display <- paste("Processing", file, sep=" ")
  print(display)
  results <- read.delim(file, header=FALSE, stringsAsFactors=FALSE)
  results <- as.numeric(as.matrix(results))
  results
  
  
  if( loop == 1) {
    results.graph <<- results

  } else {
    results.graph <<- cbind(results.graph, results)
  }
  
  loop <<- loop + 1
})

outputfilename="Unallocated_pools_size.pdf"
pdf(outputfilename, height=20, width=20)
boxplot(results.graph, main="Unallocated pool size by batch", ylab="Unallocated pool size", xaxt='n')
dev.off()




#########################################################
#Obs/sim and P-value changes with successive simulation results batches
library(data.table)


files <- list.files(path=".", pattern="\\.obs_sim_changes.txt", full.names=T, recursive=FALSE)
lapply(files, function(file) {
	changes <- fread(file)    #Read in the data using fread as conventional import is very slow
	changes <- t(changes)

	outputfilename <- paste(file, ".pdf", sep="")
	pdf(outputfilename, height=20, width=20)
	boxplot(changes, ylab="Fold change", xlab="Simulation batch", outline=FALSE, 
			main="Obs/sim changes with addition of simulation batches (outliers not shown)", cex.main=1)
	dev.off()
})



files <- list.files(path=".", pattern="\\.p_value_changes.txt", full.names=T, recursive=FALSE)
lapply(files, function(file) {
	changes <- fread(file)
	changes <- t(changes)

	outputfilename <- paste(file, ".pdf", sep="")
	pdf(outputfilename, height=20, width=20)
	boxplot(changes, ylab="Fold change", xlab="Simulation batch", outline=FALSE, 
			main="P-value changes with addition of simulation batches (outliers not shown)", cex.main=1)
	dev.off()
})




