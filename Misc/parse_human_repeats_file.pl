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
#A Perl script to parse the file human_repeats_info_sequence.txt
#to create a better annotated FASTA file
#############################################

use strict;
use warnings;

use Data::Dumper;

my $file = $ARGV[0];

my %feature_counter;    #%{feature} = count;

my $outfile = "$file.edited.fa";
open(IN, '<', $file) or die "Could not open '$file' : $!";
open(OUT, '>', $outfile) or die "Could not write  '$file' : $!";

my $id;

print "Processing $file\n";

while(<IN>){
	my $line = $_;
	chomp $line;
	
	if( substr($line,0,1) eq '>'){    #Header
		($id) = split(/\s/, $line);
		#$id =~ s/^>hg38_rmsk_//;
		
		$feature_counter{$id}++;
		
		$id = $id . '_' . $feature_counter{$id};
		print OUT ">$id\n";
	}else{    #Sequence
		print OUT "$line\n";
	}
}

	
close IN or die "Could not close $file : $!";
close OUT or die "Could not close $outfile : $!";


print "Processing complete.\n";

exit (0);

