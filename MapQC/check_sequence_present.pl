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
use Math::Round;
use Data::Dumper;
use String::Approx 'amatch';
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

#Reads a Fastq file and writes RNA reads / DNA read to sepetate files
#RNA reads have sequence GGGAGCGTGGTT starting at position 48 (1-based reference)

unless(@ARGV){
	die "Please specify a file to process.\n";
}

foreach my $file (@ARGV){

	print "Processing $file\n";

	my %counter = (Exact => 0, Absent => 0, Processed => 0, Fuzzy => 0);

    if ( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read file \'$file\' : $!";
    } elsif ( $file =~ /\.bz2$/ ) {
        open( IN, "bzcat $file |" ) or die "Couldn't read file \'$file\' : $!";
    } else {
        open( IN, $file ) or die "Couldn't read file \'$file\' : $!";
    }

	open( PRESENT, "| gzip -c - > $file.present.fastq.gz" ) or die "Couldn't write to file '$file.present.fastq.gz': $!";
	open( ABSENT, "| gzip -c - > $file.absent.fastq.gz" ) or die "Couldn't write to file '$file.absent.fastq.gz': $!";

	while (<IN>) {
			
		my $line1 = $_;
		my $line2 = scalar <IN>;
		my $line3 = scalar <IN>;
		my $line4 = scalar <IN>;
		
		$counter{Processed}++;
		
		#if($line2 =~ /^[AGCTN]{47}GGGAGCGTGGTT/){
		my $fragment = substr($line2, 26, 20);

		if($fragment =~ /GACACGCAGGGATGAGATGG/){    #Very short, use the one below in future
		#if($fragment =~ /GGGAGCGTGG/){
			print PRESENT $line1 . $line2 . $line3 . $line4;
			$counter{Exact}++;
		}elsif( amatch('GACACGCAGGGATGAGATGG', $fragment) ){
			print PRESENT $line1 . $line2 . $line3 . $line4;
			$counter{Fuzzy}++;
		}else{
			print ABSENT $line1 . $line2 . $line3 . $line4;
			$counter{Absent}++;
		}
	}

	close IN;
	close PRESENT or die "Could not close filehandle on '$file.present.fastq.gz' : $!";
	close ABSENT or die "Could not close filehandle on '$file.absent.fastq.gz' : $!";
	
	$counter{Total_Matching} = $counter{Exact} + $counter{Fuzzy};
	$counter{Percent_Matching} = calc_perc($counter{Total_Matching}, $counter{Processed}, 1);
		
	#Print summary results
	open(SUMMARY, '>', "$file.summary_sequence_present.txt") or die "Could not open '$file.summary_sequence_present.txt' : $!";
	print SUMMARY "File\tProcessed\tExact_Match\tFuzzy_Match\tNo_Match\tTotal_Matching\tPercent_Matching\n";	
	print SUMMARY "$file\t$counter{Processed}\t$counter{Exact}\t$counter{Fuzzy}\t$counter{Absent}\t$counter{Total_Matching}\t$counter{Percent_Matching}\n";
	close SUMMARY or die "Could not close filehandle to file '$file.summary_sequence_present.txt' : $!";
}

print "Processing complete.\n";

exit (0);



##########################################################
#Subroutines
##########################################################


#Subroutine: calc_perc
#Receives a number and a total and returns the perentage value
#Optional argument: decimal places of the output
# sub calc_perc {
#     my ( $n, $total, $dp ) = @_;
	
# 	if(defined $dp){
# 		$dp = abs( int($dp) );
# 	}else{
# 		$dp = 2;
# 	}
# 	$dp = 1 / (10**$dp);   #'Nearest' function needs 1, 10, 100 etc. 
	
# 	return 'NA' unless(defined $n and defined $total);   #Avoid initialisation error
# 	return 'NA' if($total == 0);    #Avoid division by zero error
	
# 	my $pc = 100 * $n / $total;
# 	$pc = nearest($dp, $pc);
	
# 	return $pc;
# }
