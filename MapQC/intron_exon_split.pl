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


##############################################################################
#Perl script takes a data file of
#Extra classification of genes:
#Nascent - any overlap with unambiguous intron
#Ejunc - both ends of split on exons
#ex= all - nasc

##############################################################################

use strict;
use warnings;
use List::Util qw( min max );
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;
use Getopt::Long;
use Data::Dumper;
use POSIX;


#Check input ok
#Option variables
my %config = (introns => undef);


my $config_result = GetOptions(
					"introns=s" => \$config{introns},
			      );

die "Could not parse options.\n" unless ($config_result);

die "Please specify files to process.\n" unless (@ARGV);
die "Please specify an --introns file to process (containing a list of unambiguous introns).\n" unless(defined $config{introns});

my @files = deduplicate_array(@ARGV);
die "Unambiguous intron file '$config{introns}' does not exist, please adjust configuration.\n" unless (-e $config{introns});
die "Please adjust configuration.\n" unless (check_files_exist(\@files, 'EXISTS'));



#Process unambiguous intron file and create intron lookup data structure
	
##############################################
#Read in capture file

print "Reading pre-specified regions file '$config{introns}'\n";


if( $config{introns} =~ /\.gz$/ ) {
    open( INTRON_FILE, "zcat $config{introns} |" ) or die "Couldn't read '$config{introns}' : $!";
} else {
    open( INTRON_FILE, '<', $config{introns} ) or die "Could not read '$config{introns}' : $!";
}


my %introns;    #%{strand}->{strand}->{chromosome}->{10kb region}->{Start\tEnd} = @[Feature ID \t Feature Name]   Expect array to have one entry

scalar <INTRON_FILE>;   #Skip header
while(<INTRON_FILE>){
	my $line = $_;
	chomp $line;
	next if $line =~ /^\s*$/;    #Ignore blank lines
	my ($csome, $start, $end, $strand) = split(/\t/, $line);
	add_bait_with_strand($csome, $start, $end, \%introns, $strand);
}
close INTRON_FILE or die "Could not close '$config{introns}' : $!";


#print Dumper \%introns;

#Create summary file
my $summary_outfile = 'intron_exon_split_summary_file.txt';
open (SUMMARY, '>', $summary_outfile) or die "Could not open '$summary_outfile' : $!";
#print SUMMARY "FILE\tREADS_PROCESSED\tEXON_JUNCTION\tNASCENT\tPUTATIVE_EXON\tOTHER\n";
print SUMMARY "FILE\tREADS_PROCESSED\tEXON_JUNCTION\tNASCENT\tOTHER\n";

#Process each file in turn
print "Processing files:\n";
foreach my $file (@files){
	print "\t$file\n";

	#Create output files
	my ($ejunc_outfile, $nascent_outfile, $all_outfile, $other_outfile) = get_filenames($file);		
	open(EXON_JUNCTION, "| samtools view -bSh - > $ejunc_outfile"  ) or die "Could not write to '$ejunc_outfile' : $!";
	open(NASCENT, "| samtools view -bSh - > $nascent_outfile"  ) or die "Could not write to '$nascent_outfile' : $!";
	#open(PUTATIVE_EXON, "| samtools view -bSh - > $all_outfile"  ) or die "Could not write to '$all_outfile' : $!";
	open(OTHER, "| samtools view -bSh - > $other_outfile"  ) or die "Could not write to '$other_outfile' : $!";
	my %summary_counter = (EXON_JUNCTION => 0, NASCENT => 0, PUTATIVE_EXON => 0, OTHER => 0);

	#Read input
	if($file =~ /\.bam$/){
		open( IN, "samtools view -h $file |" ) or die "Couldn't read '$file' : $!";
	} else {
		open (IN, '<', $file) or die "Couldn't read '$file' : $!";
	}
	
	my $header_flag = 1;
	my $header = '';
	
	while(<IN>){
		my $line = $_;
		chomp $line;

		if($header_flag){	
			if((substr($line,0,1) eq '@')){
				$header = $header . $line . "\n";
				next;
			} else {
				print EXON_JUNCTION $header;
				print NASCENT $header;
				#print PUTATIVE_EXON $header;
				print OTHER $header;
				$header_flag = 0;    #No longer in header		
			}
		}
		
		
		$summary_counter{READS_PROCESSED}++;
		
		#Definition Ejunc: Read split once, and both ends of reads map to an exon
		my (undef, $samFlag, $csome, $start, undef, $cigar, undef, undef, undef, $seq) = split(/\t/, $line);
		my $strand =  $samFlag & 0x10 ? '+' : '-';    #This is actually the wrong way round, but reads will be in an opposite orientation from "their" feature
		my $exon_junction_flag = 0;
		my $unambiguous_intron_flag = 0;
		
		#print "$samFlag\t$csome\t$start\t$cigar\t$strand\n";
		
			
		if($cigar =~ /^(\d+)M(\d+)N(\d+)M$/){    #Read spans an intron
			#Check that both ends are found in an intron (i.e. not an unambiguous intron)
			my $pos_lookup1 =  $start + int($1 / 2);   #Is this an exon (use middle to avoid imperfect splicing)
			my $pos_lookup2 = $start + $1 + $2 + int($3 / 2);
			
			if( !(coord2bin_with_strand($csome, $pos_lookup1, \%introns, $strand) ) and !(coord2bin_with_strand($csome, $pos_lookup2, \%introns, $strand)) ){
				$exon_junction_flag = 1;		
			}
		}
		
		#Definition Unambiguous intron: read not split and is entirely surrounded by an intron (which should not contain exon sequence, however the sequence is spliced)
		unless($cigar =~ /N/){ 
			my $end = $start + length($seq) - 1;
			if( (coord2bin_with_strand($csome, $start, \%introns, $strand) ) and (coord2bin_with_strand($csome, $end, \%introns, $strand)) ){
				$unambiguous_intron_flag = 1;		
			}		
		}
		
	
		#Write out the read depending on the classification
		if($exon_junction_flag and $unambiguous_intron_flag){    #Internal code check
			warn "Exon Juction Flag and Unambiguous Intron Flag both 1 at line:\n";
			warn "$line\n";
			die "This should not happen!\n";    	
		} elsif($exon_junction_flag){
			print EXON_JUNCTION "$line\n";
			#print PUTATIVE_EXON "$line\n";   #File of all putative exons	
			$summary_counter{EXON_JUNCTION}++;
			#$summary_counter{PUTATIVE_EXON}++;
		} elsif($unambiguous_intron_flag){
			print NASCENT "$line\n";
			$summary_counter{NASCENT}++;
		} else {
			#print PUTATIVE_EXON "$line\n";
			#$summary_counter{PUTATIVE_EXON}++;
			print OTHER "$line\n";
			$summary_counter{OTHER}++;
		}
	}
	close IN or die "Could not close filehandle on '$file' : $!";	
	close EXON_JUNCTION or die "Could not close filehandle on '$ejunc_outfile' : $!";
	close NASCENT or die "Could not close filehandle on '$nascent_outfile' : $!";
	#close PUTATIVE_EXON or die "Could not close filehandle on '$all_outfile' : $!";
	close OTHER or die "Could not close filehandle on '$other_outfile' : $!";

	#Write out to summary file
	#print SUMMARY "$file\t$summary_counter{READS_PROCESSED}\t$summary_counter{EXON_JUNCTION}\t$summary_counter{NASCENT}\t$summary_counter{PUTATIVE_EXON}\t$summary_counter{OTHER}\n";	
	print SUMMARY "$file\t$summary_counter{READS_PROCESSED}\t$summary_counter{EXON_JUNCTION}\t$summary_counter{NASCENT}\t$summary_counter{OTHER}\n";	
}
	
close SUMMARY or die "Could not close filehandle on '$summary_outfile' : $!";
print "Finished Exon/Intron splitting.\n";



exit (0);









###################################################################################
#Subroutines
###################################################################################


#Subroutine get_filenames
#Takes the input filename and returns the filenames of the three output files
sub get_filenames{
	my $file = $_[0];
	$file =~ s/\.filtered_reads\.bam//;
	my $ejunc_outfile = "$file.EXON_JUNCTION.filtered_reads.bam";
	my $intron_outfile = "$file.NASCENT.filtered_reads.bam";
	my $all_outfile = "$file.PUTATIVE_EXON.filtered_reads.bam";
	my $other_outfile = "$file.OTHER.filtered_reads.bam";
	return($ejunc_outfile, $intron_outfile, $all_outfile, $other_outfile); 
}






##########################################
#Subroutine: add_bait_with_strand
#Takes the bait chromosome/start/end
#and populates the passed hash accordingly:
#%{chromosome}->{10kb region}->{Start} = End
#Note: if the bin/fragment spans more than one 10kb region,
#then multiple 10 regions will be populated.
#The strand may be passed as an optional parameter
sub add_bait_with_strand {
    my ($csome, $start, $end, $hash_ref, $strand) = @_;
    	
    my $ten_kb_start = ceil($start / 10_000);
    my $ten_kb_end = ceil($end/ 10_000);
		
	die "Strand must be '+' or '-'.\n" unless($strand eq '+' or '-');	
	for (my $ten_kb_region = $ten_kb_start; $ten_kb_region <= $ten_kb_end; $ten_kb_region++){
		${$hash_ref}{$strand}->{$csome}->{$ten_kb_region}->{$start} = $end;
	}		
}



##########################################
#Subroutine: coord2bin_with_strand{
#Receives a chromosome name and a position and reference to the baits hash
#and returns the bait co-ordinates where this location is found (else returns 0)
#%lookup_hash{chromosome}->{10kb region}->{Start} = End
sub coord2bin_with_strand{
    my ($csome, $pos, $hash_ref, $strand) = @_;
    my $ten_kb_region = ceil($pos / 10_000);
	
	die "Strand must be '+' or '-'.\n" unless($strand eq '+' or '-');
	foreach my $start ( keys %{ ${$hash_ref}{$strand}->{$csome}->{$ten_kb_region} }  ){
		my $end = ${ $hash_ref }{$strand}->{$csome}->{$ten_kb_region}->{$start};
		if ( ($start <= $pos) and ($end >= $pos) ){
			return ("$csome\t$start\t$end\t$strand");
		}
	}		

    return 0;    #Not captured
}






