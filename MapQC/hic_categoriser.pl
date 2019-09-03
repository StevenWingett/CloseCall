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
#Barcode_Group	Chromosome	Position	Strand
#and returns the a a frequency table of  Hi-C category, sub-divided by
#barcode group_size
#
#Hi-C Categories
#Cis_very_close: <=1kb
#Cis_close: >1kb <=10kb   
#Cis_far: >10kb   <=10mb
#Cis_very_far: >10mb
#Trans
##############################################################################

use strict;
use warnings;
use List::Util qw( min max );
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;


#Check input ok
die "Please specify files to process.\n" unless (@ARGV);
die "Please adjust configuration.\n" unless( check_files_exist(\@ARGV, 'EXISTS') );

foreach my $file (@ARGV){

	print "Processing $file\n";
	open (IN, '<', $file) or die "Could not read '$file' : $!";

	my $outfile = "$file.hic_category_v_bargroup.txt";
	open (CATEGORY, '>', $outfile ) or die "Could not write to '$outfile' : $!";
	print CATEGORY "Group_Size\tCategory\n";
	
	#Create the random human dataset
	print "Generating random dataset\n";
	print CATEGORY "Random\t1_Cis_Wobble\n" x 1;   #Too high, but need at least 1 value
	print CATEGORY "Random\t2_Cis_Extremely_Close\n" x 1;   #Too high, but need at least 1 value
	print CATEGORY "Random\t3_Cis_Very_Close\n" x 1;   #Too high, but need at least 1 value
	print CATEGORY "Random\t4_Cis_Close\n" x 1;   #Too high, but need at least 1 value
	print CATEGORY "Random\t5_Cis_Far\n" x 6;
	print CATEGORY "Random\t6_Cis_Very_Far\n" x 43;
	print CATEGORY "Random\t7_Trans\n" x 950;
	
	
	print "Processing $file\n";
	my %groupedReads;    #%{group} = @(reads)	
	while(<IN>){
		my $line = $_;
		chomp $line;
		next if($line =~ /^\s*$/);    #Ignore whitespaces (should not be present)
		
		my ($group, $csome, $pos, $strand) = split(/\t/, $line);
		push ( @{ $groupedReads{$group} }, "$csome\t$pos" );
	}
		
	#print Dumper \%groupedReads;
	

	
	#Create the virtual HiC dataset
	foreach my $group (keys %groupedReads){
		my $group_size = scalar @{ $groupedReads{$group} };
	
		next if( scalar @{ $groupedReads{$group} } < 2 );    #Ignore n=1 barcode groups
		my @ditags = createDitags( @{ $groupedReads{$group} } );
		#print Dumper \@ditags;
		my @hic_categories = get_categories(@ditags);
		foreach my $category (@hic_categories){
			print CATEGORY "$group_size\t$category\n";
		}
	}
	
	close IN or die "Could not close '$file' : $!";
	close CATEGORY or die "Could not close '$outfile' : $!";
}

print "Processing complete.\n";

exit (0);





#######################################################################################
#Subroutines                                                                          #
#######################################################################################

#get_cis_distances
#Takes an array of ditags and returns an array of Hi-C categories
#Hi-C Categories
#Cis_very_close: <=1kb
#Cis_close: >1kb <=10kb   
#Cis_far: >10kb   <=10mb
#Cis_very_far: >10mb
#Trans
sub get_categories{
	my @ditags = @_;
	my @hic_categories;
	
	for(my $i = 0; $i < (scalar @ditags); $i+=2){
		my ($csomeA, $posA) = split(/\t/, $ditags[$i]);    #Read pairs A and B
		my ($csomeB, $posB) = split(/\t/, $ditags[$i+1]);
		my $distance = abs ($posA - $posB);

		if($csomeA ne $csomeB){
			push (@hic_categories, '7_Trans') ;    #Ignore trans
		}elsif($distance <= 10){
			push (@hic_categories, '1_Cis_Wobble');
		}elsif($distance <= 100){
			push (@hic_categories, '2_Cis_Extremely_Close');
		}elsif($distance <= 1_000){
			push (@hic_categories, '3_Cis_Very_Close');
		}elsif($distance <= 10_000){
			push (@hic_categories, '4_Cis_Close');
		}elsif($distance <= 10_000_000){
			push (@hic_categories, '5_Cis_Far');
		}else{
			push(@hic_categories, '6_Cis_Very_Far');
		}	
	}
	return @hic_categories;
}




#createDitags
#Takes an array of reads and creates an array of the read-pair combinations
#(reports read pairs once i.e. read1-read2 but not read2-read1 as well)
sub createDitags{
	my @reads = @_;
	my @ditags;
 
	while(scalar @reads > 1){
		for (my $i=1; $i < scalar @reads; $i++){
			push(@ditags, $reads[0]);
			push(@ditags, $reads[$i]);
		}	
		shift @reads;
	}
	#print Dumper \@ditags;
	return @ditags;
}




#############################################################
#check_files_exist:
#Takes a reference to an array containing paths to filenames and verifies they exist
#Warns of files that do no exit. Returns 1 if all files exist but 0 if this is not
#the case.
#
#Also, takes a second argument:
#$_[1] should be 'EXISTS' or 'NOT_EXISTS'
#If 'NOT_EXIST' warns if file already exists.  Returns '1' if none of the
#files exists and '0' if one or multiple files already exist
# sub check_files_exist {
#     my $files      = $_[0];    #Reference to array
#     my $check_for  = $_[1];
#     my $all_exist  = 1;
#     my $not_exists = 1;

#     if ( $check_for eq 'EXISTS' ) {
#         foreach my $file (@$files) {
#             unless ( -e $file ) {
#                 warn "$file does not exist\n";
#                 $all_exist = 0;
#             }
#         }
#     } elsif ( $check_for eq 'NOT_EXISTS' ) {
#         foreach my $file (@$files) {
#             if ( -e $file ) {
#                 warn "$file already exists\n";
#                 $not_exists = 0;
#             }
#         }
#     } else {
#         die "Subroutine 'check_files_exist' requires argument 'EXISTS' or 'NOT_EXISTS'.\n";
#     }

#     if ( $check_for eq 'EXISTS' ) {
#         return $all_exist;
#     } else {
#         return $not_exists;
#     }
# }

