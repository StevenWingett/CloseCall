#!/usr/bin/perl

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



