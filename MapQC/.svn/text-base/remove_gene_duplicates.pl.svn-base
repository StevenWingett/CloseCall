#!/usr/bin/perl

################################################
#Takes a file which lists features and barcode IDs where reads have and
#returns the same data but with duplicates removed (retaining 1 copy)
#
#
#
#Input Data format columns (tab-delimited):
#Read_ID 
#Barcode_ID 
#Chromosome
#Feature_Start
#Feature_End
#Feature_Strand
#Feature_ID
#Feature_Name
#
#
#Output Data format columns (tab-delimited):



################################################

use strict;
use warnings;
use Math::Round;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;

##############################################
#Check input
unless(@ARGV){
	die "Please specify an Anaconda mapped to genes file to process.\n";
}

foreach my $file (@ARGV){
	print "Processing $file\n";

	my %ids;    #%{Barcode_id\tChromosome\tStart\tEnd} = ''
	my %summary_counter = (Processed => 0, Retained => 0);
			
	#Open input file and write to output files
	if( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }

    my $outfile = $file;
    $outfile =~ s/\.on_target\.txt\.gz//;
    $outfile .= ".uniques.txt.gz";

    open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
    print OUT "Read_ID\tBarcode_ID\tChromosome\tFeature_Start\tFeature_End\tFeature_Strand\tFeature_ID\tFeature_Name\n";
    my $summary_file = "$file.uniques.summary.txt";
    open(SUMMARY, '>', "$summary_file") or die "Could not write to '$summary_file' : $!";

    scalar <IN>;    #Ignore header
    while(<IN>){
    	my $barcode_gene = $_;
    	chomp $barcode_gene;
        next if($barcode_gene =~ /^\s*$/);    #Ignore blank lines
    	$summary_counter{Processed}++;

        my @barcode_gene_id_elements = split(/\t/, $barcode_gene);
        shift @barcode_gene_id_elements;
        my $barcode_gene_id = join("\t", @barcode_gene_id_elements);

    	unless( exists $ids{$barcode_gene_id} ) {
    		print OUT "$barcode_gene\n";
    		$summary_counter{Retained}++;
    		$ids{$barcode_gene_id} = '';
    	}
    }
    close IN or die "Could not close filehandle on '$file' : $!";
    close OUT or die "Could not close filehandle on '$outfile' : $!"; 


    my $percent_retained = calc_perc($summary_counter{Retained}, $summary_counter{Processed}, 2);
    print SUMMARY "File\tProcessed\tRetained\tPercent_Retained\n"; 
    print SUMMARY "$file\t$summary_counter{Processed}\t$summary_counter{Retained}\t$percent_retained\n"; 
    close SUMMARY or die  "Could not close filehandle on '$summary_file' : $!";
}


print "Removed duplicates\n";

exit (0);

