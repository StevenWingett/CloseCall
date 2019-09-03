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


################################################################
#A Perl script that takes an Anaconda file, and reports the
#frequency of each feature by valency
#
#Anaconda file input format:
#Read_ID	
#barcode_id Chromosome	
#Feature_Start	
#Feature_End	
#Feature_Strand	
#Feature_ID	
#Feature_Name
#
#Output matrix:
#
#              Valency=1      Valency=2     Valency=3
#FeatureA         50             18             3
#FeatureB          4              0             0
#FeatureC      12034            389            51
#FeatureD         19              0             1
#
#
################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;



###########################################################################
#Verify user input

#Option variables
my %config = (
	help => ''
);
	
my $config_result = GetOptions(
    "help"     	=> \$config{help},
);
die "Could not parse options\n" unless ($config_result);

if ( $config{help} ) {    #Print help and exit
    print while (<DATA>);
    exit(0);
}

unless(@ARGV){
	warn "Please specify file(s) to process.\n";
	print while (<DATA>);
    die "\n";
}

my @files = deduplicate_array(@ARGV);


##########################################################################
#Read in data
my %features_valencies;    # %{features}->{valency} = count
my @features_in_complex;    #@ = (feature_ids}
my $max_valency = 0;
	
foreach my $file (@files){

	print "Reading in file '$file'\n";
	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	}else{
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}	
	scalar <IN>;    #Ignore header


	my $previous_barcode_id = 0;

	while(<IN>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		shift @line_elements;
		my $barcode_id = shift @line_elements;
		my $feature_name = $line_elements[5];
		
		if($barcode_id == $previous_barcode_id){    #Add feature
			push(@features_in_complex, $feature_name);	
		}else{    #Record features present in this valency
			analyse_complex();
			push(@features_in_complex, $feature_name);	
		}	
		$previous_barcode_id = $barcode_id;	
	}
	
	analyse_complex();    #May have complex at end of file
	
	#Print out results;
	my $outfile = "$file.valency_feature_distribution.txt.gz";
	open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
	print OUT "Feature";
	for(my $valency = 1; $valency <= $max_valency; $valency++){
		print OUT "\tVal=$valency";
	}
	print OUT "\n";
	
	foreach my $feature (keys %features_valencies){
		print OUT "$feature\t";
		for(my $valency = 1; $valency <= $max_valency; $valency++){
			my $frequency;
			if(exists $features_valencies{$feature}->{$valency}){
				$frequency = $features_valencies{$feature}->{$valency};
			}else{
				$frequency = 0;		
			}
			print OUT "\t$frequency";
		}
			print OUT "\n";
	}
	close OUT or die "Could not close filehandle on '$outfile' : $!";
	
	%features_valencies = ();    #File processed, empty data structure
	$max_valency = 0;    #Reset for next file
}

print "Valency distribution determined\n";

exit (0);



############################################################################################
#Subroutines
############################################################################################

# Determines valency and populates %valency_features as required
# variables not passed to subroutine, but declared outside
sub analyse_complex {
	my $valency = scalar (@features_in_complex);
	if($valency > $max_valency){
		$max_valency = $valency;
	}
	
	return unless ($valency);    #When at start and possibly end of file	
	
	foreach my $feature (@features_in_complex){
		$features_valencies{$feature}->{$valency}++;
	}
	@features_in_complex = ();    #Empty array
}













