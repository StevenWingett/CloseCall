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

#Reads a Fastq file and writes RNA reads / DNA read to sepetate files
#RNA reads have sequence GGGAGCGTGGTT starting at position 48 (1-based reference)

unless(@ARGV){
	die "Please specify a file to process.\n";
}



foreach my $file (@ARGV){

	print "Processing $file\n";

	my %counter = (DNA => 0, RNA => 0);

    if ( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read file \'$file\' : $!";
    } elsif ( $file =~ /\.bz2$/ ) {
        open( IN, "bzcat $file |" ) or die "Couldn't read file \'$file\' : $!";
    } else {
        open( IN, $file ) or die "Couldn't read file \'$file\' : $!";
    }

	open( DNA, "| gzip -c - > $file.DNA.fastq.gz" ) or die "Couldn't write to file '$file.DNA.fastq.gz': $!";
	open( RNA, "| gzip -c - > $file.RNA.fastq.gz" ) or die "Couldn't write to file '$file.RNA.fastq.gz': $!";

	open(SUMMARY, '>', "summary.$file.RNA_reads.txt") or die "Could not open 'summary.$file.RNA_reads.txt' : $!";
	print SUMMARY "File\tDNA_reads\tRNA_reads\n";	
	
	while (<IN>) {
		my $line1 = $_;
		my $line2 = scalar <IN>;
		my $line3 = scalar <IN>;
		my $line4 = scalar <IN>;
		
		#if($line2 =~ /^[AGCTN]{47}GGGAGCGTGGTT/){
		my $fragment = substr($line2, 39, 30);

		if($fragment =~ /AGCGTGG/){    #Very short, use the one below in future
		#if($fragment =~ /GGGAGCGTGG/){
			print RNA $line1 . $line2 . $line3 . $line4;
			$counter{RNA}++;
		}else{
			print DNA $line1 . $line2 . $line3 . $line4;
			$counter{DNA}++;
		}
	}

	close IN;
	close DNA or die "Could not close filehandle on '$file.DNA.fastq.gz' : $!";
	close RNA or die "Could not close filehandle on '$file.RNA.fastq.gz' : $!";

	print SUMMARY "$file\t$counter{DNA}\t$counter{RNA}\n";

}

close SUMMARY or die "Could not close filehandle to file 'summary_RNA_reads.txt' : $!";

print "Processing complete.\n";

exit (0);