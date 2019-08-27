#!/usr/bin/perl

#############################################
#A Perl script to parse the file human_repeats_info_sequence.txt
#to create a better annotated FASTA file
#############################################

use strict;
use warnings;

use Data::Dumper;

my $outfile = "ncRNA.edited.fa";
open(OUT, '>', $outfile) or die "Could not write  '$outfile' : $!";

foreach my $file (@ARGV){

	if ($file =~ /^Homo_sapiens\.GRCh38\.(.+)\.fa$/){
		my $id = $1;
		my $i = 1;    #%{feature} = count;

		print "Processing $file\n";
		open(IN, '<', $file) or die "Could not open '$file' : $!";
		while(<IN>){
			my $line = $_;
			chomp $line;
			
			if( substr($line,0,1) eq '>'){    #Header
				print OUT ">$id" . '_' . "$i\n";
				$i++;
			}else{    #Sequence
				print OUT "$line\n";
			}
		}		
		close IN or die "Could not close $file : $!";
	}else{
		warn "Skipping '$file' : does not have correct file name structure: /^Homo_sapiens\.GRCh38\.(.+)\.fa$/\n";
	}
}
close OUT or die "Could not close $outfile : $!";

print "Processing complete.\n";

exit (0);

