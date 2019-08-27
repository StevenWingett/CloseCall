#!/usr/bin/perl

#############################################
#A Perl script to combine Anaconda Simulation format
#files.  This is achieved by adding 10,000,000 to the
#Read ID and Barcode ID of subsequent files.
#i.e.:
#FILE1
#Read_ID Barcode_ID      Chromosome      Feature_Start   Feature_End     Feature_Strand  Feature_ID      Feature_Name
#1       1       1       11908151        11908271        +       353     U5
#2       2       11      62853733        62854160        -       54188   32766.SNHG1
#
#
#FILE2
#Read_ID Barcode_ID      Chromosome      Feature_Start   Feature_End     Feature_Strand  Feature_ID      Feature_Name
#4       4       14      20323179        20323326        -       58826   39015.SNORA79
#5       5       gi|374429547|ref|NR_046235.1|   7925    12994   +       46602   6.28S
#
#becomes:
#COMBINED FILE
#1       1       1       11908151        11908271        +       353     U5
#2       2       11      62853733        62854160        -       54188   32766.SNHG1
#100000004       100000004       14      20323179        20323326        -       58826   39015.SNORA79
#100000005       100000005       gi|374429547|ref|NR_046235.1|   7925    12994   +       46602   6.28S
#############################################

use strict;
use warnings;

use Data::Dumper;

unless(@ARGV){
        die "Please specify files to process.\n";
}

my @files = deduplicate_array(@ARGV);
my $outfile = 'anaconda_simulation_format_combined_file.txt.gz';
open(OUT, "| gzip -c - > $outfile" ) or die "Couldn't write to file '$outfile' : $!";
my $i = 0;
print "Processing:\n";
foreach my $file (@files){

	print "\t$file\n";

	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	} else {
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}
	
	my $header = scalar <IN>;
	print OUT $header unless($i);

	while(<IN>){
		my $line = $_;
		chomp $line;
		my @data = split(/\t/, $line);
		$data[0] = $data[0] + ($i * 100_000_000);
		$data[1] = $data[1] + ($i * 100_000_000);
		print OUT join("\t", @data);
		print OUT "\n";
	}
	
	$i++;
	close IN or die "Could not close '$file' : $!";
}
close OUT or die "Could not close '$outfile' : $!";
  
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







