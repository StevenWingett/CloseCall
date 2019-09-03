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


use strict;
use warnings;
use Data::Dumper;


#Reads a Fastq file and fixes the barcode to 26 As

unless(@ARGV){
	die "Please specify a file to process.\n";
}

foreach my $file (@ARGV){

	print "Processing $file\n";

    if ( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read file \'$file\' : $!";
    } elsif ( $file =~ /\.bz2$/ ) {
        open( IN, "bzcat $file |" ) or die "Couldn't read file \'$file\' : $!";
    } else {
        open( IN, $file ) or die "Couldn't read file \'$file\' : $!";
    }

	open( OUT, "| gzip -c - > $file.fixed.fastq.gz" ) or die "Couldn't write to file '$file.fixed.fastq.gz': $!";
	

	while (<IN>) {
			
		my $line1 = $_;
		my $line2 = scalar <IN>;
		my $line3 = scalar <IN>;
		my $line4 = scalar <IN>;

		$line2 = substr($line2, 26);
		my $polyA = 'A' x 26;
		$line2 = $polyA . $line2;
		print OUT $line1 . $line2 . $line3 . $line4;
	}

	close IN;
	close OUT or die "Could not close filehandle on '$file.fixed.fastq.gz' : $!";
	
}

print "Processing complete.\n";

exit (0);

