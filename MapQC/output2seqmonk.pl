#!/usr/bin/perl

use strict;
use warnings;

###############################################################################################
#Script to convert file containing Hi-C data in a format produced by the Monte Carlo script:
#Read1 Csome, Read1 Start, Read1 End,  Read2 Csome, Read2 Start, Read2 End   Count
#
#into a format compatible with SeqMonk i.e.
#Read1 Csome    Read1 Start    Read1 End
#Read2 Csome    Read2 Start    Read2 End
#.
#With count resulting in the same read represented multiple times
##############################################################################################

unless (@ARGV) {
    die "Please enter a filename";
}

foreach my $file (@ARGV) { 
    print "Processing $file\n";

    if ( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
    } else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }

	my $outfile = $file . '_seqmonk.txt.gz';
	open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
	print OUT "Chromosome\tStart\tEnd\n";
    scalar <IN>;    #Ignore header line

    while (<IN>) {
	
		my $line = $_;
		chomp $line;

        my ( $chromosomeF, $startF, $endF, $chromosomeR, $startR, $endR, $count ) = split(/\t/, $line);

        #Convert Chromosome 'M' to 'MT'
        if ( $chromosomeF eq 'M' ) {
            $chromosomeF = 'MT';
        }
        if ( $chromosomeR eq 'M' ) {
            $chromosomeR = 'MT';
        }
        for ( my $i = 1 ; $i <= $count ; $i++ ) {
            print OUT "$chromosomeF\t$startF\t$endF\n$chromosomeR\t$startR\t$endR\n";
        }
    }
}

print "Processing complete\n";

exit;

