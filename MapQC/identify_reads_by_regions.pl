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


#############################################
#Takes a condensed format file from Anaconda along with a list
#of genomic regions and returns reads that map to those regions
#Return the co-ordinates of the target region, not the original
#read (when mapping to a target region)
#############################################

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Math::Round;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;


##############################################
#Option variables
my $read_length = 50;    #Length of read mapped by HISAT2

my %config = (
	baits => '',
	header => 1,
	help => '',
	interactions => ''
);
	
my $config_result = GetOptions(
	"baits=s"	=> \$config{baits},
	"header=i" => \$config{header},
    "help"     	=> \$config{help},
	"interactions" => \$config{interactions},
);
die "Could not parse options\n" unless ($config_result);

if ( $config{help} ) {    #Print help and exit
    print while (<DATA>);
    exit(0);
}

print "Identifing reads that map to pre-specified regions\n";

##############################################
#Check input
unless(@ARGV){
warn "Please specify an Anaconda condensed format file to process.\n";
    print while (<DATA>);
    die "\n";
}
die "Please specify a --baits regions file (tab-delimited: chr  start  end) (see --help for more details)\n" if $config{baits} eq '';
die "Please adjust configuration.\n" unless ( check_files_exist(\@ARGV, 'EXISTS') );
die "Capture file '$config{baits}' does not exits, please adjust configuration.\n" unless ( -e $config{baits} );


##############################################
#Read in capture file

print "Reading pre-specified regions file '$config{baits}'\n";


if( $config{baits} =~ /\.gz$/ ) {
    open( BAITS, "zcat $config{baits} |" ) or die "Couldn't read '$config{baits}' : $!";
} else {
    open( BAITS, '<', $config{baits} ) or die "Could not read '$config{baits}' : $!";
}


if($config{header}){    #Number of lines to skip
	$config{header} = abs($config{header});    #Ensure +ve
	for (1..$config{header}){
		scalar <BAITS>;
	}
}

my %baits;    #%{strand}->{chromosome}->{10kb region}->{Start\tEnd} = @[Feature ID \t Feature Name]   Expect array to have one entry
my %bait_ids;    #%{chr\tstart\tend} = id

my $feature_id = 1;
while(<BAITS>){
	my $line = $_;
	chomp $line;
	next if $line =~ /^\s*$/;    #Ignore blank lines
	my ($csome, $start, $end, $strand, $feature_name) = split(/\t/, $line);
	add_bait_features($csome, $start, $end, $strand, $feature_id, $feature_name, \%baits);
	$feature_id++;
}
close BAITS or die "Could not close '$config{baits}' : $!";


#print Dumper \%baits;

##############################################
#Read in condensed Anaconda files and identify captured/uncaptured reads
foreach my $file (@ARGV){

	print "Processing $file\n";

	my %bait_bait_interactions;    #%{bait1_chr\tbait1_start\tbait1_end\tbait2_chr\tbait2_start\tbait2_end} = count

	#Categorise and record the data
	#initialise_hash(\%bait_bait_interactions);    #Sets all terms to 0
	
	my %results_counter;
	my @results_categories = ('Reads_Processed', 'On_Target', 'Off_Target', 'Partially_Off_Target', 'Ambiguous',
		'Percent_On_Target');
	
	foreach my $category (@results_categories){
		$results_counter{$category} = 0;
	}
			
	#Open input file and write to output files
	if( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }
	

	my $on_target_file = $file;
	$on_target_file =~ s/\.condensed_format\.txt$//;
	$on_target_file .= '.on_target.txt.gz';
	my $off_target_file = $file;
	$off_target_file =~ s/\.condensed_format\.txt$//;
	$off_target_file .= '.off_target.txt.gz';

	my $summary_filename = "$file.on_target_summary.txt";
	open(ON_TARGET, "| gzip -c - > $on_target_file") or die "Could not write to '$on_target_file' : $!";
	print ON_TARGET "Read_ID\tBarcode_ID\tChromosome\tFeature_Start\tFeature_End\tFeature_Strand\tFeature_ID\tFeature_Name\n";

	open(OFF_TARGET, "| gzip -c - > $off_target_file") or die "Could not write to '$off_target_file' : $!";
	open(SUMMARY, '>', $summary_filename) or die "Could not open '$summary_filename' : $!";
	
	{    #Print header to summary file
		my $summary_header = "File\t";
		foreach my $category (@results_categories){
			$summary_header = $summary_header . "$category\t";
		}
		$summary_header =~ s/\t$/\n/;
		print SUMMARY $summary_header;
	}


	while(<IN>){
		my $read = $_;
		chomp $read; 
		my ($barcode_group, $csome, $pos, $strand) = split(/\t/, $read);
		if($strand eq '+'){   #This is an antisense library so we need to reverse the signs !!!
			$strand = '-';
		}elsif($strand eq '-'){
			$strand = '+';
		}
				
		$results_counter{Reads_Processed}++;

		my @regions = coord2bin_features($csome, $pos, $strand, \%baits);
		if (scalar @regions > 1){    #Only one regions should be returned
			warn "Read mapping to multiple regions: '$read'\n";
			warn "Read maps to regions:\n";
			foreach my $region (@regions){
				warn "$region\n";		
			}
			die "Check the features input file.  Do the repeats overlap with one another?\n";
		}

		if(@regions){    #On target
			
			#Check the start, middle and end of read all map to same feature
			my $midpoint_region = $regions[0];
		#	my @five_prime_pos_offset_regions = coord2bin_features($csome, ($pos - int($read_length/2)), $strand, \%baits);    #Only check 1 region, don't check for multiple regions here
		#	my @three_prime_pos_offset_regions = coord2bin_features($csome, ($pos + int($read_length/2)), $strand, \%baits); 
				
				#if(@five_prime_pos_offset_regions and @three_prime_pos_offset_regions){    #Non-empty arrays returned
				#	if( ($five_prime_pos_offset_regions[0] eq $midpoint_region) and ($three_prime_pos_offset_regions[0] eq $midpoint_region) ){ 
						print ON_TARGET "$results_counter{Reads_Processed}\t$barcode_group\t$midpoint_region\n";    #Print out the target region
						$results_counter{On_Target}++;
				#	}else{
				#		print OFF_TARGET "$read\n";     #Print out the read
				#		$results_counter{Ambiguous}++;    
				#	}				
				#}else{
				#	print OFF_TARGET "$read\n";     #Print out the read
				#	$results_counter{Partially_Off_Target}++;    		
				#}
		
		}else{    #Off target
			print OFF_TARGET "$read\n";     #Print out the read
			$results_counter{Off_Target}++;    
		}
	}

	$results_counter{Percent_On_Target} = calc_perc( $results_counter{On_Target}, $results_counter{Reads_Processed} );
	{   #Print results to the summary file
		my $summary_results = "$file\t";
		foreach my $category (@results_categories){
			$summary_results = $summary_results . "$results_counter{$category}\t";
		}
		$summary_results =~ s/\t$/\n/;
		print SUMMARY $summary_results;
	}
	
	close IN or warn "Could not close '$file' : $!";
	close ON_TARGET or warn "Could not close '$on_target_file' : $!";
	close OFF_TARGET or warn "Could not close '$off_target_file' : $!";

}

print "Processing complete.\n";

exit (0);




##########################################################################
#Subroutines
##########################################################################


##########################################
#Subroutine: add_bait_features
#Takes the bait chromosome/start/end
#and populates the passed hash accordingly:
#%{strand}->{chromosome}->{10kb region}->{Start\tEnd} = "Feature_ID\tFeature_Name"   !Expect to have one value
#Note: if the bin/fragment spans more than one 10kb region,
#then multiple 10 regions will be populated
sub add_bait_features {
	my ($csome, $start, $end, $strand, $feature_id, $feature_name, $hash_ref) = @_;
		
	my $ten_kb_start = ceil($start / 10_000);
	my $ten_kb_end = ceil($end/ 10_000);
	
	for (my $ten_kb_region = $ten_kb_start; $ten_kb_region <= $ten_kb_end; $ten_kb_region++){
		if(exists ${$hash_ref}{$strand}->{$csome}->{$ten_kb_region}->{"$start\t$end"}  ){
			warn "Feature $strand $csome $start $end present more than once\n";
		} 
		${$hash_ref}{$strand}->{$csome}->{$ten_kb_region}->{"$start\t$end"} = "$feature_id\t$feature_name";
	}
}




##########################################
#Subroutine: coord2bin_features
#Receives a chromosome name and a position and reference to the baits hash
#and returns the bait co-ordinates where this location is found (else returns empty array)
#%{strand}->{chromosome}->{10kb region}->{Start\tEnd} = "Feature_ID\tFeature_Name"   !Expect to have one value
sub coord2bin_features{
	my ($csome, $pos, $strand, $hash_ref) = @_;

	my $ten_kb_region = ceil($pos / 10_000);
	my @results;

	foreach my $start_end ( keys %{ ${$hash_ref}{$strand}->{$csome}->{$ten_kb_region} } ){
		my ($start, $end) = split(/\t/, $start_end);
		if ( ($start <= $pos) and ($end >= $pos) ){
			my ($feature_id, $feature_name) = split( /\t/, ${$hash_ref}{$strand}->{$csome}->{$ten_kb_region}->{$start_end} );

			push(@results, "$csome\t$start\t$end\t$strand\t$feature_id\t$feature_name");
		}
	}

	return @results;
}

#${$hash_ref}{$strand}->{$csome}->{$ten_kb_region}->{"$start\t$end"} 
__DATA__

