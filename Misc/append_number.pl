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


#############################################
#Append 1 billion to fist column in file
#############################################

use strict;
use warnings;

use Data::Dumper;

unless(@ARGV){
	die "Please specify a file to process.\n";
}

foreach my $file(@ARGV){
	my $outfile = "appended.$file";
	open(IN, '<', $file) or die "Could not open '$file' : $!";
	open(OUT, '>', $outfile) or die "Could not write to '$outfile' : $!";

	while(<IN>){
		my $line = $_;
		chomp $line;
		my @data = split(/\t/, $line);
		$data[0] = $data[0] + 1_000_000_000;
		print OUT join("\t", @data);
		print OUT "\n";
	}

	
	close IN or die "Could not close $file : $!";
	close OUT or die "Could not close $outfile : $!";


}


print "Processing complete\n";

exit (0);
