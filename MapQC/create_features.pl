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
#A Perl script to take a features file and
#parse it so overlapping genes ON THE SAME STRAND become 
#separate entities.  Regions of overlap are NOT assigned
#a feature UNLESS the feature is contained completely within
#another feature.
#A strand may either be +/-/B, where B means it is to considered
#as both sense and antisense

#1         1   1000      +      GeneA
#1      2000   2100      +      GeneB
#1      2050   2500      +      GeneC
#1         1   3000      -      GeneD
#2         1   1000      +      GeneE
#2       101    200      +      GeneF
#3          1   700      B       Rep1      
#becomes:
#
#1         1   1000      +      GeneA     1.1
#1      2000   2049      +      GeneB     2.2
#1      2101   2500      +      GeneC     3.3
#1         1   3000      -      GeneD     4.4
#2         1    100      +      GeneE     5.5
#2       101    200      +      GeneF     6.6
#2       201   1000      +      GeneE     5.7
#3         1    700      +       Rep1     7.8
#3         1    700      -       Rep1     7.9

#############################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;
use POSIX;

use Data::Dumper;

###################################################################
#Input options
my %config = (
	header => 1,
	help => '',
);
	
my $config_result = GetOptions(
	"header=i" => \$config{header},
    "help"     	=> \$config{help},
);
die "Could not parse options\n" unless ($config_result);

if ( $config{help} ) {    #Print help and exit
   # print while (<DATA>);
    exit(0);
}

unless(@ARGV){
warn "Please specify a feature file to process.\n";
    #print while (<DATA>);
    die "\n";
}


#####################################################
#Read input file and write data to variables
my $file= shift @ARGV;

if( $file =~ /\.gz$/ ) {
    open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
} else {
    open( IN, '<', $file ) or die "Could not read ' $file' : $!";
}


if($config{header}){    #Number of lines to skip
	$config{header} = abs($config{header});    #Ensure +ve
	for (1..$config{header}){
		scalar <IN>;
	}
}

#Read input features file and store in a relevant data structure
#my %points_of_change;   #%{strand}->{chromosome}->{position} = ''
my %features;    #%{$strand}->{$csome}->{Start/End?}->{$pos}-> = @(ID\Name)
my %repeats;    # %{"+\t$csome\t$start\t$end\t$id\t$name"}

print "Processing $file\n";

my $id = 1;
my $repeats_encountered_flag = 0;    #A check to help ensure input file in correct format
while(<IN>){
	my $line = $_;
	chomp $line;
	next if($line =~ /^\s*$/);

	my ($csome, $start, $end, $strand, $name) = split(/\t/, $line);
	$end = $end + 1;       #This is the actual point where the annotation changes
	
	if($strand eq 'B'){    #Patch: add twice to both strands to ensure this region is free from features
		push( @{ $features{'+'}->{$csome}->{$start}->{Start} }, "$id\t$name");
		push( @{ $features{'+'}->{$csome}->{$end}->{End} }, "$id\t$name" );
		push( @{ $features{'-'}->{$csome}->{$start}->{Start} }, "$id\t$name");
		push( @{ $features{'-'}->{$csome}->{$end}->{End} }, "$id\t$name" );
		
		push( @{ $features{'+'}->{$csome}->{$start}->{Start} }, "$id\t$name");
		push( @{ $features{'+'}->{$csome}->{$end}->{End} }, "$id\t$name" );
		push( @{ $features{'-'}->{$csome}->{$start}->{Start} }, "$id\t$name");
		push( @{ $features{'-'}->{$csome}->{$end}->{End} }, "$id\t$name" );
		
		$repeats{"+\t$csome\t$start\t$end\t$id\t$name"} = '';
		$repeats{"-\t$csome\t$start\t$end\t$id\t$name"} = '';

		$repeats_encountered_flag = 1;
			
	}else{
		if($repeats_encountered_flag){    #Problem in input file
			warn "The input file '$file' should contain repeats (denoted by strand is B) only AFTER all other features\n";
			die "Terminating.\n";
		}
		push( @{ $features{$strand}->{$csome}->{$start}->{Start} }, "$id\t$name");
		push( @{ $features{$strand}->{$csome}->{$end}->{End} }, "$id\t$name" );
	}
	$id++;
	
}
close IN or die "Could not close filehandle on '$file' : $!";

%repeats = merge_overlapping_repeats(%repeats);    #Resolve problems in input file - let's make sure the input file is correct instead

my $results_string = '';

#Now loop through the data structure 'features' and record all the non-overlapping structures 
#%{$strand}->{$csome}->{$pos}->{Start/End}-> = @(ID\Name)
my @id_names_in_play;    #Stores ids started by not ended (most recent feature first)
my $looking_for_start = 1;   #Flag determining whether to print start
foreach my $strand (keys %features){
	foreach my $csome ( keys %{ $features{$strand} } ){
		foreach my $pos (  sort {$a <=> $b} keys %{ $features{$strand}->{$csome} } ){  

			####End - remove a feature###
			if(exists $features{$strand}->{$csome}->{$pos}->{End}){
				my %id_names_to_remove;    #%{id} = name
				my @updated_id_names_in_play;
				foreach my $id_name ( @{ $features{$strand}->{$csome}->{$pos}->{End} } ){
					$id_names_to_remove{$id_name} = '';
				}

				foreach my $id_name (@id_names_in_play){
					push(@updated_id_names_in_play, $id_name) unless (exists $id_names_to_remove{$id_name});
				}

				#End feature if only 1 feature in play
				if( (scalar @id_names_in_play == 1) and (scalar keys %id_names_to_remove) ){    
					my $actual_end = $pos - 1;
					$results_string .= "$actual_end\t$id_names_in_play[0]\n";
				}

				#Ending a feature may cause a new feature to start
				if( (scalar @updated_id_names_in_play == 1) and (scalar @id_names_in_play != 1) ){
					$results_string .= "$strand\t$csome\t$pos\t"
				}

				@id_names_in_play = @updated_id_names_in_play;
			}



			####Start - add to array###
			if(exists $features{$strand}->{$csome}->{$pos}->{Start}){
				my @id_names_to_add;
				foreach my $id_name ( @{ $features{$strand}->{$csome}->{$pos}->{Start} } ){    #Get the new ids to add		
					push (@id_names_to_add, $id_name);
				}


				#Start new feature with nothing in play
				if(scalar @id_names_to_add == 1){    #Only 1 new feature starts here
					$results_string .= "$strand\t$csome\t$pos\t" unless(@id_names_in_play);  #Already in a feature  						
				}

				#If 1 feature already in play, end previous feature
				if(scalar @id_names_in_play == 1){
					my $end_previous_feature = $pos - 1;
					$results_string .= "$end_previous_feature\t$id_names_in_play[0]\n";    
				}

				push(@id_names_in_play, @id_names_to_add);
			}
		}
	}
}



#Now loop again through the data structure 'features' and record all 
#fragments that start and end with no other fragments starting or ending
#in between.
foreach my $strand (keys %features){
	foreach my $csome ( keys %{ $features{$strand} } ){
		my $current_id_name = 'Not_Valid';
		my $current_start_pos;
		my $current_end_pos;
		foreach my $pos (  sort {$a <=> $b} keys %{ $features{$strand}->{$csome} } ){  
			
			#Start of a feature
			if(exists $features{$strand}->{$csome}->{$pos}->{Start}){
				my @id_names;
				$current_start_pos = $pos;
				foreach my $id_name ( @{ $features{$strand}->{$csome}->{$pos}->{Start} } ){    #Get the new ids to add		
					push (@id_names, $id_name);
				}

				if(scalar @id_names == 1){    #Only 1 feature start at this point
					$current_id_name = $id_names[0];
				}else{
					$current_id_name = 'Not_Valid';
					next;
				}

				#Check no feature ends at this point
				if(exists $features{$strand}->{$csome}->{$pos}->{End}){
					$current_id_name = 'Not_Valid';
					next;
				}
			}


			#End of a feature
			if(exists $features{$strand}->{$csome}->{$pos}->{End}){
				my @id_names;
				my $current_end_id_name;
				$current_end_pos = $pos - 1;
				foreach my $id_name ( @{ $features{$strand}->{$csome}->{$pos}->{End} } ){    #Get the new ids to add		
					push (@id_names, $id_name);
				}

				if(scalar @id_names == 1){    #Only 1 feature ends at this point
					$current_end_id_name = $id_names[0];   #_END_ID_NAME
				}else{
					$current_id_name = 'Not_Valid';
					next;
				}

				#Check no feature start at this point
				if(exists $features{$strand}->{$csome}->{$pos}->{Start}){
					$current_id_name = 'Not_Valid';
					next;
				}

				#Is this an uninterrupted feature
				if( ($current_id_name eq $current_end_id_name) and ($current_id_name ne 'Not_Valid') ){
					$results_string .= "$strand\t$csome\t$current_start_pos\t$current_end_pos\t$current_id_name\n";
				}else{
					$current_id_name = 'Not_Valid';
				}

				$current_id_name = $current_id_name;

			}
		}
	}
}

#print Dumper \@id_names_in_play;
#print "$results_string\n\n\n";

#Remove duplicate values from the results string and print results
my $outfile = basename($file);
$outfile .= ".edited_features.txt.gz";

open(OUT, "| gzip -c - > $outfile") or die "Couldn't write to file '$outfile' : $!";
print OUT "Chromosome\tStart\tEnd\tStrand\tFeature_ID\n";
my @results = split(/\n/, $results_string);
my %unique_results;


#Remove any of these features found entirely within a repeat region 
#(Feature partially overlapping should already have been removed)
#Create a lookup hash of the repeat regions
# %{chromosome}->{10kb region}->{Start\tEnd} = ''
my %repeats_lookup_hash;
foreach my $repeat (keys %repeats){
	my (undef, $csome, $start, $end) = split(/\t/, $repeat);    #Use the key, not the hash value
	add_repeat_bait($csome, $start, $end, \%repeats_lookup_hash);
}

foreach my $result (@results){
	my (undef, $csome, $start, $end) = split(/\t/, $result);
	unless(coord2bin_repeat($csome, $start, $end, \%repeats_lookup_hash) ) {  #Found WITHIN a repeat
		$unique_results{$result} = '';	
	}	
}


#Create a separate file of the modified repeats
# %{"+\t$csome\t$start\t$end\t$id\t$name"}
my $modified_repeats_file = basename($file);
$modified_repeats_file .= ".edited_repeats_only.txt.gz";
open(MODIFIED_REPEATS, "| gzip -c - > $modified_repeats_file") or die "Couldn't write to file '$modified_repeats_file' : $!";
print MODIFIED_REPEATS "Chromosome\tStart\tEnd\tStrand\tRepeat_ID\n";
foreach my $repeat (keys %repeats){
	my($strand, $csome, $start, $end, $id, $name) = split(/\t/, $repeat);
	$end = $end - 1;   #End correction adjustment!!!
	print MODIFIED_REPEATS "$csome\t$start\t$end\t$strand\t$name\n";
}
close MODIFIED_REPEATS or die "Could not close filehandle on '$modified_repeats_file' : $!";
print "Created repeats file '$modified_repeats_file'\n";


foreach my $result (keys %repeats){    #Add repeats
	#Remove 1 from the end of the repeat
	my ($strand, $csome, $start, $end, $id, $name) = split(/\t/, $result);
	$end = $end - 1;    #End correction
	my $new_result = join("\t", ( $strand, $csome, $start, $end, $id, $name ) );
	$unique_results{$new_result} = '';
}

my %feature_id_counter;

#print Dumper \%unique_results;




foreach my $result (sort keys %unique_results){   #Print results and append count index to feature id
	my @feature_id_elements = split(/\t/, $result);
	my $feature_id = $feature_id_elements[4];
	$feature_id_counter{$feature_id}++;
	my $count = $feature_id_counter{$feature_id};
	$feature_id_elements[4] = "$feature_id_elements[4].$count";


	#Rearrange outfile, so it resembles BED format and merge ID components
	my $feature_name = pop(@feature_id_elements);
	$feature_id_elements[-1] = "$feature_id_elements[-1].$feature_name";
	my @part2 = splice(@feature_id_elements, 0, 1);
	my @part1 = splice(@feature_id_elements, 0, 3);

	push(@part1, @part2, @feature_id_elements);
	my $result = join("\t", @part1);
	print OUT "$result\n";
}

close OUT or die "Could not close filehandle on '$outfile' : $!";

print "Created repeats/features file '$outfile'\n";

exit (0);






##########################################################################
#Subroutines
##########################################################################

###########################################
#Subroutine: merge_overlapping_repeats
#The input file may not be perfect with regard to repeats i.e. repeats may overlap.
#To resolve this merge overlapping repeats (which map to the same repeat class)/
#Takes the repeat hash an returns a modified repeats hash.
#%repeats{Strand Chromosome    Start    End    Id    Name} = '';
sub merge_overlapping_repeats{
	my %repeats = @_;
	#print Dumper \%repeats;
	
	
	my %repeats_restructured;   #%{Chromosome}->{Start End} = {Id    Name}  #Strand can be ignored
	
	foreach my $repeat (keys %repeats){    #Create new data structure
		my (undef, $csome, $start, $end, $id, $name) = split(/\t/, $repeat);
		
		if(exists $repeats_restructured{$csome}->{$start}->{$end}){
			unless($repeats_restructured{$csome}->{$start}){    #Check for problem
				warn "Different repeat features start at same position....edit input file:\n";
				die "Chromsome:$csome Start: $start\n";	
			}
		}
		$repeats_restructured{$csome}->{$start}->{$end} = "$id\t$name";
	}
	
	#print Dumper \%repeats_restructured;
	
	%repeats = ();    # Clear hash
	
	foreach my $csome (keys %repeats_restructured){
		my $previous_start = 0;
		my $previous_end = 0;
		my $previous_id;
		my $previous_name;
		my $i = 0;
		foreach my $start ( sort {$a <=> $b} keys %{ $repeats_restructured{$csome} } ){
			foreach my $end ( sort {$a <=> $b} keys %{ $repeats_restructured{$csome}->{$start} } ){    #Should only be one value
				$i++;	
				my $id_name = $repeats_restructured{$csome}->{$start}->{$end};
				my ($id, $name) = split(/\t/, $id_name);
				
				#Does this overlap?
				if($start > $previous_end){    #Not overlap
					unless($i == 1) {    #First feature is chromosome
						$repeats{"+\t$csome\t$previous_start\t$previous_end\t$previous_id\t$previous_name"} = '';
						$repeats{"-\t$csome\t$previous_start\t$previous_end\t$previous_id\t$previous_name"} = '';
					}
					
					$previous_start = $start;
					$previous_end = $end;
					$previous_id = $id;
					$previous_name = $name;
						
				}else{    #Overlaps	
					if($name eq $previous_name){
						$previous_end = $end;	
					}else{
						die "Repeat name overlap conflict at chromsome $csome at $start-$end and $previous_start-$previous_end.\n";		
					}
				}	
			}
		}
	
		$repeats{"+\t$csome\t$previous_start\t$previous_end\t$previous_id\t$previous_name"} = '';
		$repeats{"-\t$csome\t$previous_start\t$previous_end\t$previous_id\t$previous_name"} = '';	

	}
	
	#print Dumper \%repeats;
	return %repeats;
	exit;
}








##########################################
#Subroutine: add_repeat_bait
#Takes the bait chromosome/start/end
#and populates the passed hash accordingly:
# %{chromosome}->{10kb region}->{Start\tEnd} = ''   
#Note: if the bin/fragment spans more than one 10kb region,
#then multiple 10 regions will be populated
#Strand is not important since repeats will be allocated to both strands
sub add_repeat_bait {
	my ($csome, $start, $end, $hash_ref) = @_;
		
	my $ten_kb_start = ceil($start / 10_000);
	my $ten_kb_end = ceil($end/ 10_000);
	
	for (my $ten_kb_region = $ten_kb_start; $ten_kb_region <= $ten_kb_end; $ten_kb_region++){
		${$hash_ref}{$csome}->{$ten_kb_region}->{"$start\t$end"} = '';
	}
}




##########################################
#Subroutine:  coord2bin_repeat
#Check wheter feature is found WITHIN a repeat
#%{chromosome}->{10kb region}->{Start\tEnd} = ''
sub coord2bin_repeat{
	my ($csome, $lookup_start, $lookup_end, $hash_ref) = @_;

	my $ten_kb_region = ceil($lookup_start / 10_000);
	my @results;

	foreach my $start_end ( keys %{ ${$hash_ref}{$csome}->{$ten_kb_region} } ){
		my ($start, $end) = split(/\t/, $start_end);
		if ( ($lookup_start >= $start) and ($lookup_end <= $end) ){
			return 1;
		}
	}
	
	return 0;
}




