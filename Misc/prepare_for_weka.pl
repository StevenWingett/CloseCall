#!/usr/bin/perl

#################################################################
#Perl script to prepare a file for Weka Random Forest#
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
#Return in Weka format the interacting features and attributes that interact with (and only with)
#the features of interest
#####################################################################


use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;

##################################
#Get user input

#Option variables
my $default_minimum_count = 3;
my $default_minimum_p_value = 0.01;


my %config;
my $config_result = GetOptions(
								"target=s" => \$config{target},
								"matrix=s" => \$config {matrix},
								"simulation=s" => \$config {simulation},
								"count=i" => \$config {count},
								"p=f" => \$config{p}		   
							);

die "Could not parse options" unless ($config_result);

unless( (defined $config{target}) and (defined $config{matrix}) and (defined $config{simulation}) ){
	die "Please specify --target, --matrix and a --simulation file.\n";
}

unless(defined $config{count}){
	warn "Setting minumum --count to $default_minimum_count\n";
	$config{count} = $default_minimum_count;
}

unless(defined $config{p}){
	warn "Setting minumum --p value to $default_minimum_p_value\n";
	$config{p} = $default_minimum_count;
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

#print Dumper \%valid_interacting_features;


################################################################
# Annotate each feature with attributes and write out the result
if($config{matrix} =~ /\.gz$/){
	open (MATRIX, "zcat $config{matrix} |") or die "Could not read file '$config{matrix}' : $!";
} else {
	open(MATRIX, '<', $config{matrix}) or die "Could not open '$config{matrix}' : $!";
}

my $outfile = 'random_forest_weka_input.csv.gz';
open(OUT, "| gzip -c - > $outfile") or die "Couldn't write to file '$outfile' : $!";
my $header = scalar <MATRIX>;
chomp $header;
$header =~ s/\t/,/g;    #Comma-separated output
print OUT "$header,Class\n";

while(<MATRIX>){
	my $line = $_;
	chomp $line;
	
	my ($interacting_feature) = split(/\t/, $line);
	if(exists $valid_interacting_features{$interacting_feature}){
		my $output = "$line\t$valid_interacting_features{$interacting_feature}\n";   #Print out line and class
		$output =~ s/\t/,/g; #Comma-separated output
		print OUT $output;
		delete $valid_interacting_features{$interacting_feature};    #Remove from hash, and check if hash empty at the end
	}
}
close MATRIX or die "Could not close '$config{matrix}' : $!";
close OUT or die "Could not close '$outfile' : $!";


if(scalar %valid_interacting_features){
	print "Warning: the following remain in the 'valid_interacting_features' hash:\n";
		foreach my $interacting_feature (keys %valid_interacting_features){
			print "\t$interacting_feature\n";
		}
}


print "Processing complete\n";
exit (0);