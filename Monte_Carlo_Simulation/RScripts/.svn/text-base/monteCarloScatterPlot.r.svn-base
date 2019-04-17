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