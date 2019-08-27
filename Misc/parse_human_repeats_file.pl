#!/usr/bin/perl

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

