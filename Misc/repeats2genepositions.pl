#!/usr/bin/perl

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


#############################################################
#A Perl script for parsing a file so
#Chromosome	Start	End	Name	Class	Family	Anaconda_Name
#1	67108753	67109046	L1P5	LINE	L1	L1
#
#is parsed to:
#
#Chromosome	Start	End	Strand	ID
#1	11869	14409	+	DDX11L1
#
#
#Every feature in the input file will have a + and - Strand
#in the output file
############################################################

use strict;
use warnings;

my ($file) = @ARGV;

if ($file =~ /\.gz$/){
	open (IN,"zcat $file |") or die "Couldn't read $file : $!";  
}else{
	open (IN, $file) or die "Could not open $file\n";
}

my $outfile .= "$file.parsed2gene_format.txt";
open(OUT, '>', $outfile) or die "Couldn't write to '$outfile' : $!";

scalar <IN>;   #Ignore header

while(<IN>){
	my $line = $_;
	chomp $line;

	my ($csome, $start, $end, undef, undef, undef, $name) = split(/\t/, $line);

	print OUT "$csome\t$start\t$end\t+\t$name\n";
	print OUT "$csome\t$start\t$end\t-\t$name\n";    #Both strands
}

close IN;
close OUT or die "Could not close filehandle out '$outfile' : $!";

print "Done\n";

exit (0); 



