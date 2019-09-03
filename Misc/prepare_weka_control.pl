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


#################################################################
#Perl script to prepare a control file for Weka Random Forest#
#Takes:
#1) A --target list of genes of interest (1 column file) e.g. 
#Neat1
#U3
#
#2) A --matrix of genes/features and associated attributes (e.g. TF present)
#Records 0/1 to denote whether an attribute is present
#Feature CCNT2   BRF2    Brg1    c-Fos   c-Jun
#GeneA
#GeneB
#GeneC
#
#3) The results of an Anaconda Monte Carlo --simulation
#
#4) You can also specify the minimum absolute interaction --count allowed [default 3]
#
#5) You can also specify the maximum allowed interaction --p value
#
#6) It is also possible to create a control file which contains the feature that DO NOT
#interact with the pre-specified list
# 
#Return in Weka format the interacting features and attributes that interact with (and only with)
#the features of interest
#####################################################################


use strict;
use warnings;
use Getopt::Long;
use List::Util qw(shuffle);

use Data::Dumper;

##################################
#Get user input

#Option variables
my $default_minimum_count = 3;
my $default_maximum_p_value = 0.05;


my %config;
my $config_result = GetOptions(
								"name=s" => \$config{name},
								"target=s" => \$config{target},
								"matrix=s" => \$config {matrix},
								"simulation=s" => \$config {simulation},
								"count=i" => \$config {count},							
								"p=f" => \$config{p}		   
							);

die "Could not parse options" unless ($config_result);

unless( (defined $config{target}) and (defined $config{matrix}) and (defined $config{simulation}) and (defined $config{name} ) ){
	die "Please specify --target, --matrix, --simulation file and a --name\n";
}

if(defined $config{count}){
	warn "Minimum --count set by user to $config{count}\n";
}else{
	$config{count} = $default_minimum_count;
	warn "Setting minumum --count to $config{count}\n";
}

if(defined $config{p}){
	warn "Maximum --p value set by user to $config{p}\n";
}else{
	$config{p} = $default_maximum_p_value;
	warn "Setting minumum --p value to $config{p} \n";

}


###############################
#Identify features of interest
my %features;    # %{feature_of_interest}->{interacting_feature} = ''
if($config{target} =~ /\.gz$/){
	open (TARGET, "zcat $config{target} |") or die "Could not read file '$config{target}' : $!";
} else {
	open(TARGET, '<', $config{target}) or die "Could not open '$config{target}' : $!";
}

while(<TARGET>){
	my $line = $_;
	chomp $line;
	$features{$line} = undef;    #Set as undef to prevent autovivification error later on
}
close TARGET or die "Could not close '$config{target}' : $!";


#######################################################
#Read simulation file to identify interacting features
if($config{simulation} =~ /\.gz$/){
	open (SIM, "zcat $config{simulation} |") or die "Could not read file '$config{simulation}' : $!";
} else {
	open(SIM, '<', $config{simulation}) or die "Could not open '$config{simulation}' : $!";
}


my %control_features;   # %control{feature} = ''

scalar <SIM>;    #Ignore header
while(<SIM>){
	my $line = $_;
	chomp $line;
	my ($feat1, $feat2, $count, undef, undef, undef, $pval) = split(/\t/, $line);
	$feat1 =~ s/^\d+\.//;    #Remove numerical prefix
	$feat2 =~ s/^\d+\.//;
	
	next if($count < $config{count});    #Check meets parameters
	next if($pval > $config{p});
	
	my $feature_of_interest;
	my $interacting_feature;
	
	if( (exists $features{$feat1}) and (exists $features{$feat2}) ){    #Ignore if an interaction between 2 features of interest
		next;
	}elsif(exists $features{$feat1}){
		$feature_of_interest = $feat1;
		$interacting_feature = $feat2;
	}elsif(exists $features{$feat2}){	
		$feature_of_interest = $feat2;
		$interacting_feature = $feat1;
	}else{    #Not of interest
		$control_features{$feat1} = '';
		$control_features{$feat2} = '';
		next;
	}	
	$features{$feature_of_interest}->{$interacting_feature} = ''; 
}

close SIM or die "Could not close '$config{simulation}' : $!";






#Identify and remove interacting features that interact with more than one feature of interest
my %valid_interacting_features;  #%{feature} = count;
foreach my $feature_of_interest (keys %features){    
	if(defined $features{$feature_of_interest}){
		foreach my $interacting_feature ( keys %{ $features{$feature_of_interest} } ){
			$valid_interacting_features{$interacting_feature}++;
		}			
	}
}

foreach my $interacting_feature (keys %valid_interacting_features){
	if($valid_interacting_features{$interacting_feature} > 1){
		delete $valid_interacting_features{$interacting_feature};   #Remove from hash
	}
}

#Obtain the original feature of interest
foreach my $feature_of_interest (keys %features){    
	if(defined $features{$feature_of_interest}){
		foreach my $interacting_feature ( keys %{ $features{$feature_of_interest} } ){
			if(exists $valid_interacting_features{$interacting_feature}){
				$valid_interacting_features{$interacting_feature} = $feature_of_interest;
			}
		}
	}
}



my %valid_interacting_features_copy = %valid_interacting_features;    #No elements to be deleted from these
my %control_features_copy = %control_features;
#print Dumper \%valid_interacting_features;


################################################################
# Annotate each feature with attributes and write out the result
if($config{matrix} =~ /\.gz$/){
	open (MATRIX, "zcat $config{matrix} |") or die "Could not read file '$config{matrix}' : $!";
} else {
	open(MATRIX, '<', $config{matrix}) or die "Could not open '$config{matrix}' : $!";
}

my $outfile = $config{name} . '.random_forest_weka_input.csv';
my $controlfile = $config{name} . '.random_forest_weka_control.csv';
my $randContDataFile = $config{name} . '.random_forest_data_random_control_subset.csv';


open(OUT, '>', $outfile) or die "Couldn't write to file '$outfile' : $!";
open(CONTROL, '>', $controlfile) or die "Couldn't write to file '$controlfile' : $!";
open(DATACON, '>', $randContDataFile) or die "Couldn't write to file '$randContDataFile' : $!";


my $header = scalar <MATRIX>;
chomp $header;
$header =~ s/\t/,/g;    #Comma-separated output
print OUT "$header,Class\n";
print CONTROL "$header,Class\n";
print DATACON "$header,Class\n";

while(<MATRIX>){
	my $line = $_;
	chomp $line;
	
	my ($interacting_feature) = split(/\t/, $line);
	if( (exists $valid_interacting_features{$interacting_feature}) ){    #Is this a control dataset or a real dataset
		my $output = "$line\t$valid_interacting_features{$interacting_feature}\n";   #Print out line and class
		$output =~ s/\t/,/g; #Comma-separated output
		print OUT $output;
		print DATACON $output;
		delete $valid_interacting_features{$interacting_feature};    #Remove from hash, and check if hash empty at the end
	} elsif (exists $control_features{$interacting_feature} ){	
		my $output = "$line\n";   #Print out line and class
		$output =~ s/\t/,/g; #Comma-separated output
		print CONTROL $output;
		delete $control_features{$interacting_feature};
		#delete $valid_interacting_features{$interacting_feature};    #Remove from hash, and check if hash empty at the end
	}
}
close MATRIX or die "Could not close '$config{matrix}' : $!";
close OUT or die "Could not close '$outfile' : $!";
close CONTROL or die "Could not close '$controlfile' : $!";



#Report features to which attributes could not be allocated
if(scalar (keys %valid_interacting_features) ){
	print "Note: Attributes could not be assigned to the features:\n";
	foreach my $interacting_feature (keys %valid_interacting_features){
		print "\t$interacting_feature\n";
	}
}

if(scalar (keys %control_features) ){
	print "Note: Attributes could not be assigned to the following control features:\n";
	foreach my $interacting_feature (keys %control_features){
		print "\t$interacting_feature\n";
	}
}


#########################################################################
#Now make a final file with the data and the control data (random subset)

foreach my $control_feature (keys %control_features){    # Remove feature that could not be allocated attribute
	if(exists $control_features{$control_feature} ){	
		delete $control_features_copy{$control_feature};
	}
}

foreach my $feature (keys %valid_interacting_features){    # Remove feature that could not be allocated attribute
	if(exists $valid_interacting_features{$feature} ){	
		delete $valid_interacting_features_copy{$feature};
	}
}


my $number_of_features_to_use = scalar (keys %valid_interacting_features_copy);
if(  $number_of_features_to_use > scalar (keys %control_features_copy) ){
	die "Fewer features in control than data!!\n";
}

my %random_control_subset;
my @random_subset = keys(%control_features_copy);
@random_subset = shuffle(@random_subset);
@random_subset = splice(@random_subset, 0, $number_of_features_to_use);

foreach my $feature (@random_subset){
	$random_control_subset{$feature} = '';
}


if($config{matrix} =~ /\.gz$/){
	open (MATRIX, "zcat $config{matrix} |") or die "Could not read file '$config{matrix}' : $!";
} else {
	open(MATRIX, '<', $config{matrix}) or die "Could not open '$config{matrix}' : $!";
}

while(<MATRIX>){
	my $line = $_;
	chomp $line;
	my ($interacting_feature) = split(/\t/, $line);
	
	if( (exists $random_control_subset{$interacting_feature}) ){    #Is this a control dataset or a real dataset
		my $output = "$line\tCONTROL\n";   #Print out line and class
		$output =~ s/\t/,/g; #Comma-separated output
		print DATACON $output;	
	}
}

close DATACON or die "Could not close '$randContDataFile' : $!";
close MATRIX or die "Could not close '$config{matrix}' : $!";

print "Processing complete\n";

exit (0);

