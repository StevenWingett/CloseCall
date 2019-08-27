#! /usr/bin/perl
use strict;
use warnings;

# this script is for seaprating the different types of ncRNA for a genome downloaded from Ensembl.
# creates separate files for each type of ncRNA.
# http://www.ensembl.org/info/data/ftp/index.html

my $input_file = shift @ARGV;

unless ($input_file =~ /.fa|.fasta$/){
	warn "Input must be a fasta file\n";
}	

open (IN,$input_file) or die "Can't read $input_file: $!";

my %filehandles;
my $sequence = "";
my $id = "";
my $transcript_biotype = "";
my $counter = 0;

while (my $line = <IN>){

	chomp $line;
	$counter ++;
	
	if($line =~ /^>/){

		# print out the last sequence to the file handle
		if($counter > 1){
			
			unless (defined $filehandles{$transcript_biotype}) {
				
				my $outfile = $input_file;
				$outfile =~ s/\.fa|\.fasta/_$transcript_biotype.fa/;
				warn "Writing to $outfile\n";
				open (my $fh,"> $outfile") or die "Can't write to $outfile: $!";
			    $filehandles{$transcript_biotype} = $fh;
			}
		
			
			print {$filehandles{$transcript_biotype}} "$id\n";
			print {$filehandles{$transcript_biotype}} "$sequence\n"; 
			$sequence="";
		}
	
		$id = $line;
		my @split_id = split(/ /,$line);
		my $last_field = $split_id[-1];
		my @last_fields = split(/:/,$last_field);
		unless($last_fields[0] =~ /transcript_biotype/){ 
			warn "Couldn't find transcript biotype for $line \n";
		}	
		$transcript_biotype = $last_fields[1];
	}	
		
	# the sequence is on separate lines.
	else {
		$sequence = $sequence.$line;
	}
}
# write out last line
print {$filehandles{$transcript_biotype}} "$id\n";
print {$filehandles{$transcript_biotype}} "$sequence\n"; 

close IN;
