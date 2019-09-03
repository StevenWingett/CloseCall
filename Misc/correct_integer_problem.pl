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

my $nos_sims = 10_000;
my $overflows = 2;
my $int_max = 2147483647;
my $int_min = -2147483648;



unless(@ARGV){
		die "Please specify files to process.\n";
}
my @files = deduplicate_array(@ARGV);


my $correction_value = ($int_max + abs($int_min)) * $overflows;
print "Processing:\n";
foreach my $file (@files){

	print "\t$file\n";

	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	} else {
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}
	
	my $outfile = "$file.edited.txt.gz";
	open(OUT, "| gzip -c - > $outfile") or die "Couldn't write to file '$outfile' : $!";
	
	while(<IN>){
		my $line = $_;
		chomp $line;
		
		my @elements = split(/\t/, $line);
		my $feature1 = $elements[0];
		my $feature2 = $elements[1];
		my $observed = $elements[2];
		my $sim_freq = $elements[3];
		my $obs_sim = $elements[4];
		
		if( ($feature1 eq '2.18S') and ($feature2 eq '6.28S') ){
			my $sim_total = $sim_freq * $nos_sims;
			my $true_total = $sim_total + $correction_value;
			my $true_sim_average = $true_total / $nos_sims;
			
			$elements[3] =  $true_sim_average;   #Set new sim score;
			$elements[4] = $observed / $true_sim_average;
			
			print OUT join("\t", @elements) . "\n";	
		} else {
			print OUT "$line\n";
		}
	}
	close IN or die "Could not close '$file' : $!";
	close OUT or die "Could not close '$outfile' : $!";
  
}


print "Processing complete\n";

exit (0);




######################################################################
#Sub: deduplicate_array
#Takes and array and returns the array with duplicates removed
#(keeping 1 copy of each unique entry).
sub deduplicate_array{
		my @array = @_;
		my %uniques;

		foreach (@array){
				$uniques{$_} = '';
		}
		my @uniques_array = keys %uniques;
		return @uniques_array;
}
