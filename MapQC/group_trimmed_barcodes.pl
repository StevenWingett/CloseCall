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


##############################################################
#Takes a barcode SAM file (made by map_barcodes.pl
#and produces a file relating each barcode to a barcode group
#Barcode \t Group
##############################################################

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;



#Option variables
my $virtual_genome_filename;
my $config_result = GetOptions(
						"virtual=s" => \$virtual_genome_filename,
			      );

die "Could not parse options" unless ($config_result);

unless(defined $virtual_genome_filename){
	die "Please specify a --virtual genome *.fa file.\n";

}


unless(@ARGV){
	die "Please specify a file to process.\n";
}


#Read in the virtual genome
print "Reading virtual genome file: '$virtual_genome_filename'\n";

open(VIRTUAL, '<', $virtual_genome_filename) or die "Could not open '$virtual_genome_filename' : $!";
scalar <VIRTUAL>;    #Ignore header

my $barcode_count = 0;
my $barcode_length = 0;
my $max_barcodes_per_element = 1_000_000_000 / 50;    #1 billion bp / 50 bps
my @virtual_genome;   #Virtual barcode genome chunked into 1 billion bps (get 64-bit Perl problems if too large)
while(<VIRTUAL>){
	my $line = $_;
	chomp $line;

	$barcode_count++;
	my $chunk = ceil($barcode_count / $max_barcodes_per_element);
	$virtual_genome[$chunk - 1] .= $line;

}	
close VIRTUAL or die "Could not close '$virtual_genome_filename' : $!";





#Create summary file
my @summary_categories = qw(File Sequences_Processed Contains_Indels Reverse_Complement_Alignment Contains_N
							Sequences_Grouped Unique_Barcodes Unique_Barcode_Groups );

foreach my $file (@ARGV){

	open(SUMMARY, '>', "$file.group_barcodes_summary.txt") or die "Could not write to '$file.group_barcodes_summary.txt' : $!";
	my $summary_header = '';
	foreach my $category ( @summary_categories  ){
		$summary_header .= "$category\t";
	}
	$summary_header =~ s/\t$/\n/;    #Replace last tab with new line
	print SUMMARY $summary_header;


	#Create and initialise summary counter
	my %summary_counter;
		foreach my $category (@summary_categories){
		$summary_counter{$category} = 0;
	}
	$summary_counter{File} = $file;
	

	print "Processing $file\n";
	open(IN, '<', $file) or die "Could not open '$file' : $!";
	
	my %group_barcode;    #%{barcode_group}->{all barcodes in same group}
	my %barcode_group;    #%{barcode}->{group to which barcode belongs}
	
	while(<IN>){
	
		my $line = $_;
		chomp $line;
		next if(substr($line, 0, 1) eq '@');
		
		$summary_counter{Sequences_Processed}++;
				
		my @line_elements = split(/\t/, $line);
		my ($index, $barcode) = split(/\./, $line_elements[0]);
		my $sam_flag = $line_elements[1];
		my $cigar =  $line_elements[5];
		my $conversion_string = $line_elements[-2];   #?????
		
	
		unless($cigar =~ /^20M$/){  #Ignore alignments with insertions/deletions
			$summary_counter{Contains_Indels}++;
			next;
		}
		
		if( $conversion_string =~ /\^/ ){
			$summary_counter{Contains_Indels}++;
			next;
		}

		if($sam_flag & 0x10){  #Ignore reverse complement alignments
			$summary_counter{Reverse_Complement_Alignment}++;
			next;
		}
	

		#my $barcode_genome_sequence = get_genomic_sequence($barcode, $conversion_string);
		my $barcode_genome_sequence = get_genomic_sequence($line, \@virtual_genome);
		
		#print "Line: $line\n";
		#print "Barcode: $barcode_genome_sequence\n\n";
			
		if($barcode_genome_sequence =~ /N/){   #Ignore alignments containing an N
			$summary_counter{Contains_N}++;
			next;
		}
				

		$summary_counter{Sequences_Grouped}++;
		
		#Populate barcode_hash
		barcode_grouper( \%group_barcode, \%barcode_group, $barcode, $barcode_genome_sequence );
				
	}
	
	close IN or die "Could not close '$file'\n";

	#Write results to file
	open(RESULTS, '>', "$file.barcode_groups.txt") or die "Could not write to  '$file.barcode_groups.txt' : $!";
	foreach my $barcode	(keys %barcode_group){

		my $group = $barcode_group{$barcode};
		print RESULTS "$barcode\t$group\n";
	}	
	close RESULTS or die "Could not close '$file.barcode_groups.txt' : $!";
	

	#Produce summary graphs
	my %freq_counter;  
	my %large_groups;   #%{group} = size;
	foreach my $group (keys %group_barcode){
		my $group_size = scalar keys %{ $group_barcode{$group} };
		$freq_counter{$group_size}++;	
		$large_groups{$group} = $group_size if($group_size > 20);
	}

	#Large groups written to file
	open(LARGE_GROUPS, '>', "$file.barcode_large_groups.txt") or die "Could not write to  '$file.barcode_large_groups.txt' : $!";
	foreach my $value (sort {$large_groups{$a} <=> $large_groups{$b} } keys %large_groups) {
     	print LARGE_GROUPS "$value\t$large_groups{$value}\n";
	}
	close LARGE_GROUPS or die "Could not close '$file.barcode_large_groups.txt' : $!";

	open(FREQ, '>', "$file.barcode_groups_frequency.txt") or die "Could not write to  '$file.barcode_groups_frequency.txt' : $!";
	foreach my $group_size (sort {$a <=> $b} (keys %freq_counter)){
		print FREQ "$group_size\t$freq_counter{$group_size}\n";
	}	
	close FREQ or die "Could not close '$file.barcode_groups_frequency.txt' : $!";
	
	
	#Print results to summary file
	$summary_counter{Unique_Barcodes} = scalar keys %barcode_group;
	$summary_counter{Unique_Barcode_Groups} = scalar keys %group_barcode;

	my $summary_line = '';
	foreach my $category (@summary_categories ){
		$summary_line .= "$summary_counter{$category}\t";
	}
	$summary_line =~ s/\t$/\n/;    #Replace last tab with new line
	print SUMMARY $summary_line;
	
}

close SUMMARY or die "Could not write to 'group_barcodes_summary.txt' : $!";

print "Processing complete.\n";


exit (0);


###################################################################
#Subroutines
###################################################################


#Subroutine: barcode_grouper
#Takes ref to %{barcode_group}->{all barcodes in same group}, ref to %{barcode}->{group to which barcode belongs}
#and barcode and barcode_genome_sequence and updates the 2 hashes accordingly
sub barcode_grouper{
	my ( $group_barcode_ref, $barcode_group_ref, $barcode, $barcode_genome_sequence ) = @_;
		
	if( exists ${ $barcode_group_ref }{$barcode} and exists ${ $barcode_group_ref }{$barcode_genome_sequence} ){    #Both belong to group
		
		my $barcode_group_to_join = ${ $barcode_group_ref }{$barcode};
		my $barcode_group_to_delete = ${ $barcode_group_ref }{$barcode_genome_sequence};

		return if($barcode_group_to_join eq $barcode_group_to_delete);   #Barcodes already belong to same group  

		#Merge barcode groups	
		my @barcodes_to_change_group = keys %{${ $group_barcode_ref }{$barcode_group_to_delete} };
		delete ${ $group_barcode_ref }{$barcode_group_to_delete} ;   #Remove group

		 foreach my $barcode_change (@barcodes_to_change_group){
			 ${ $barcode_group_ref }{$barcode_change} = $barcode_group_to_join;    #Changed group
			 $ { $group_barcode_ref }{$barcode_group_to_join}->{$barcode_change} = '';	#Add to new group
		 }

	}elsif (exists ${ $barcode_group_ref }{$barcode} ){   #One belongs to a group
	
		#Orphan barcode joins same group as partner
		my $barcode_group_to_join = ${ $barcode_group_ref }{$barcode};
		${ $barcode_group_ref }{$barcode_genome_sequence} = $barcode_group_to_join;     #%{barcode} = group 
		$ { $group_barcode_ref }{$barcode_group_to_join}->{$barcode_genome_sequence} = '';    # %{group}->{barcode} = '';

	}elsif( exists ${ $barcode_group_ref }{$barcode_genome_sequence}  ) {   #The other belongs to a group
	
		#Orphan barcode joins same group as partner
		my $barcode_group_to_join = ${ $barcode_group_ref }{$barcode_genome_sequence};
		${ $barcode_group_ref }{$barcode} = $barcode_group_to_join;     #%{barcode} = group 
		$ { $group_barcode_ref }{$barcode_group_to_join}->{$barcode} = '';    # %{group}->{barcode} = '';
	
	}else{    #Neither belongs to a group
	
		#Create a new group
		my $barcode_group_to_join = $barcode;
		${ $barcode_group_ref }{$barcode} = $barcode_group_to_join;    #%{barcode} = group
		${ $barcode_group_ref }{$barcode_genome_sequence} = $barcode_group_to_join;
		
		$ { $group_barcode_ref }{$barcode_group_to_join}->{$barcode} = '';    # %{group}->{barcode} = '';
		$ { $group_barcode_ref }{$barcode_group_to_join}->{$barcode_genome_sequence} = ''; 
	}

}


#Subroutine: get_genomic_sequence
#Takes a SAM alignment (of the 20bp trimmed barcode) and the virtual genome sequence and 
#returns the 26bp barcode sequence to where the read mapped
sub get_genomic_sequence{
	my ($alignment, $virtual_genome_ref) = @_;
	my (undef, $midpoint) = samMidPoint($alignment);
	#print "midpoint: $midpoint\n";
	
	my $chunk = ceil($midpoint / 1_000_000_000);
	my $midpoint_pos_on_chunk = $midpoint - (1_000_000_000 * ($chunk - 1) );    #Correct position for chunk
	my $barcode_nos = ceil($midpoint_pos_on_chunk / 50);


	my $barcode = substr( ${$virtual_genome_ref}[$chunk - 1], (($barcode_nos - 1) * 50), 26 );    #Get the barcode sequence
	return $barcode;
}
