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
