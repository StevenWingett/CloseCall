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
#separate entities

#1         1   1000      +      GeneA
#1      2000   2100      +      GeneB
#1      2050   2500      +      GeneC
#1         1   3000      -      GeneD
#
#becomes:
#
#1         1   1000      +      GeneA
#1      2000   2049      +      GeneB
#1      2050   2100      +      GeneB/GeneC
#1      2101   2500      +      GeneC
#1         1   3000      -      GeneD   
#############################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

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
my %features_add_remove;    #%{csome}->{$start}->{$end}->{$strand}->{Add/remove} = @{features}

print "Processing $file\n";

while(<IN>){
	my $line = $_;
	chomp $line;
	my ($csome, $start, $end, $strand, $id) = split(/\t/, $line);
	
	$end = $end + 1;    #This is the actual point where the annotation changes
	
	push( @{ $features_add_remove{$strand}->{$csome}->{$start}->{'add'} }, $id );
	push( @{ $features_add_remove{$strand}->{$csome}->{$end}->{'remove'} }, $id );
	
}
close IN or die "Could not close filehandle on '$file' : $!";

#print Dumper \%features_add_remove;

my @results_in_blocks;
foreach my $strand ( sort { $a cmp $b } keys %features_add_remove){
	foreach my $csome ( sort { $a cmp $b } keys %{ $features_add_remove{$strand} } ) {
		my %features_in_play;
		my @results;
		foreach my $pos ( sort {$a <=> $b} keys %{ $features_add_remove{$strand}->{$csome} } ){
			
			if(exists $features_add_remove{$strand}->{$csome}->{$pos}->{add}  ){    #Add first			
				#Add new values
				my @new_features = @{ $features_add_remove{$strand}->{$csome}->{$pos}->{add} };			
				incrememntHash(\%features_in_play, 'add', \@new_features);

			}

			if(exists $features_add_remove{$strand}->{$csome}->{$pos}->{remove}  ){    #Add first			
				#Remove values
				my @remove_features = @{ $features_add_remove{$strand}->{$csome}->{$pos}->{remove} };			
				incrememntHash(\%features_in_play, 'remove', \@remove_features);
			}
			
			my @features = keys (%features_in_play);
			my $feature_string = convert2FeatureString(@features);
			
			push(@results, "$pos\t$feature_string");
			
			
		}
			push(@results_in_blocks, convert2blocks(\@results, $csome, $strand) );
					#print Dumper \%features_in_play;
	}	

}



my $outfile = basename("$file.edited_features.txt.gz");
open(OUT, "| gzip -c - > $outfile") or die "Couldn't write to file '$outfile' : $!";
print OUT "Chromosome\tStart\tEnd\tStrand\tFeature\tFeature_ID\n";
my $i = 1;
foreach my $result (@results_in_blocks){
	print OUT "$result\tFeature" . $i . "\n";
	$i++;
}
close OUT or die "Could not close '$outfile' : $!";


#########################################################
#Subroutines
#########################################################

sub convert2blocks{
	my ($input_array_ref, $csome, $strand) = @_; 
	
	my $array_size = scalar @{$input_array_ref};
	
	my @output;
	for (my $i = 0; $i < $array_size - 1 ; $i++ ){   #Loop to penultimate value
		my ($start, $feature) = split( /\t/, ${ $input_array_ref }[$i] );
		my ($end) = split(/\t/, ${$input_array_ref}[$i+1]);
		$end = $end - 1;
		push(@output, "$csome\t$start\t$end\t$strand\t$feature") unless ($feature eq '');		
	}

	return @output;
}


sub convert2FeatureString {
	my @features = deduplicate_sort_array(@_);
	my $feature_string = join('/', @features);
	return $feature_string;
}


sub deduplicate_sort_array{
	my @array = @_;
	my %uniques;

	foreach (@array){
		$uniques{$_} = '';	
	}
	my @uniques_array = sort { $a cmp $b } keys %uniques;
	return @uniques_array;
}	



sub incrememntHash {
	my ($hash_ref, $command, $features_ref) = @_;
	
	if($command eq 'add'){
		foreach my $feature ( @{ $features_ref } ){
			${ $hash_ref }{$feature}++;
		}
	}else{
		foreach my $feature ( @{ $features_ref } ){
			${ $hash_ref }{$feature}--;
		}	
		foreach my $feature ( @{ $features_ref } ){
			delete ${ $hash_ref }{$feature} if ((exists ${ $hash_ref }{$feature} ) and (${ $hash_ref }{$feature} == 0) );
		}
	}
}







