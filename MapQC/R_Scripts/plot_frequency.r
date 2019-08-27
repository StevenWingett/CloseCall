args <- commandArgs(TRUE)
file <- args[1]
title_prefix <- args[2]
outfile_suffix <- args[3]

outputfilename=paste(file, outfile_suffix, sep = ".")
pdf(file=outputfilename)

tally <- read.delim(file, header=TRUE)
graphTitle <- paste( title_prefix, file,  sep = ": ")
y_max <- max(tally$Tally)
y_max <- ceiling(y_max * 1.25)
plot( tally, type="l", ylim=c(1,y_max), main=graphTitle, cex.main=0.75 )

dev.off()
