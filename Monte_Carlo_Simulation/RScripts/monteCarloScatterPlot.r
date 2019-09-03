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


#Wait for data to be written to file
Sys.sleep(120)


#Plot Scatter plot of simulation 

files <- list.files(path=".", pattern="\\.qval.txt.gz", full.names=T, recursive=FALSE)

lapply(files, function(file) {
  
  
  file <- "test.sim.txt.simulation_collated_data.txt.gz.qval.txt.gz"
  results <- read.delim(file, header=TRUE, stringsAsFactors = FALSE)
  results <- subset(results, results$Observed_Frequency > 1)   #Only show features observed more than twice
  
  #For graphical purposes and to avoid -Inf errors, swap 0 for 0.0000001
  results$Simulation_Average_Frequency <- replace(results$Simulation_Average_Frequency, 
                                                   results$Simulation_Average_Frequency == 0, 0.0000001)
  
  
  xmax <- log((max(results$Simulation_Average_Frequency) * 1.2), 10)
  ymax <- log((max(results$Observed_Frequency) * 1.2), 10)
  
  graph_res <- cbind(results$Simulation_Average_Frequency, results$Observed_Frequency, results$p.normal.)
  
  title <- paste("Interacting Features Count Scatter Plot (observed vs simulated):", file, "Showing observed > 1", 
                 "Blue: Monte Carlo P-value < 0.001", "Green: Pass Benjamini Hochberg", sep="\n")
  

  outputfilename=paste(file, "pdf", sep = ".")
  pdf(outputfilename, height=20, width=20)
  
  plot( log(graph_res[,1:2], 10),
        xlim=c(0, xmax), ylim=c(0, ymax),   xaxs="i", yaxs="i", xlab="Log10(simulated)", ylab="Log10(observed)", 
        main=title, cex.lab=1.5, cex.axis=1.5)
  
  abline(0,1)
  cex.main=3
  
  
  #Colour by P-value
  par(new = TRUE)
  
  graph_res_pval <- subset(graph_res, graph_res[,3] <= 0.001)
  plot( log(graph_res_pval[,1:2], 10), col="blue", xlim=c(0, xmax), ylim=c(0, ymax),   
        xaxs="i", yaxs="i", xlab="Log10(simulated)", ylab="Log10(observed)", 
        main=title, cex.lab=1.5, cex.axis=1.5)
  
  
  
  
  #Colour by Benjamini-Hochberg threshold
  par(new = TRUE)
  graph_res_qval <- subset(results, results$Passed_threshold == "pass", select=c(Simulation_Average_Frequency, Observed_Frequency))
  
  plot( log(graph_res_qval[,1:2], 10), col="green", xlim=c(0, xmax), ylim=c(0, ymax),
        xaxs="i", yaxs="i", xlab="Log10(simulated)", ylab="Log10(observed)", 
        main=title, cex.lab=1.5, cex.axis=1.5)
  
  dev.off()

})