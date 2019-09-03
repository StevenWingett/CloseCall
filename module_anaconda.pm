package module_anaconda;
require Exporter;

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


our @ISA    = qw (Exporter);
our @EXPORT = qw(VERSION hasval deduplicate_array check_files_exist calc_perc
                samMidPoint add_bait coord2bin deduplicate_array);
our @EXPORT_OK = qw();


our $VERSION = "0.1.0";


use strict;
use warnings;
use Math::Round;
use POSIX;
#use File::Basename;
#use FindBin '$Bin';
#use lib $Bin;

use Data::Dumper;

##################################################################
#A collection of Perl subroutines for the Anaconda MapQC Pipeline#
##################################################################

##########################################
#Subroutine: check_files_exist
#Takes a reference to an array containing paths to filenames and verifies they exist
#Warns of files that do no exit. Returns 1 if all files exist but 0 if this is not
#the case.
#
#Also, takes a second argument:
#$_[1] should be 'EXISTS' or 'NOT_EXISTS'
#If 'NOT_EXISTS' warns if file already exists.  Returns '1' if none of the
#files exists and '0' if one or multiple files already exist
sub check_files_exist {
    my $files      = $_[0];    #Reference to array
    my $check_for  = $_[1];
    my $all_exist  = 1;
    my $not_exists = 1;

    if ( $check_for eq 'EXISTS' ) {
        foreach my $file (@$files) {
            unless ( -e $file ) {
                warn "File '$file' does not exist\n";
                $all_exist = 0;
            }
        }
    } elsif ( $check_for eq 'NOT_EXISTS' ) {
        foreach my $file (@$files) {
            if ( -e $file ) {
                warn "File '$file' already exists\n";
                $not_exists = 0;
            }
        }
    } else {
        die "Subroutine 'check_files_exist' requires argument 'EXISTS' or 'NOT_EXISTS'.\n";
    }

    if ( $check_for eq 'EXISTS' ) {
        return $all_exist;
    } else {
        return $not_exists;
    }
}




###############################################################
#Sub: hasval
#Takes a string and returns true (i.e. '1') if the string has a value
#(i.e. is not equal to nothing (''). This is useful since some
#variables may be set to nothing allowing them to be evaluated
#without initialisation errors, while simultaneously containing
#no information.
sub hasval {
    if ( $_[0] ne '' ) {
        return 1;
    } else {
        return 0;
    }
}



##########################################
#Subroutine: samMidPoint
#Receives a SAM format line and returns
#chromosome, midpoint of the position of then
#read
sub samMidPoint {
    my $read = $_[0];

    my $csome = ( split( /\t/, $read ) )[2];
    my $start_genome_perspective = ( split( /\t/, $read ) )[3];
    my $seq   = ( split( /\t/, $read ) )[9];
    
    my $length = length($seq);  
    my $midpoint = $start_genome_perspective + ceil($length / 2);

    return ( $csome, $midpoint );
}


##########################################
#Subroutine: add_bait
#Takes the bait chromosome/start/end
#and populates the passed hash accordingly:
#%{chromosome}->{10kb region}->{Start} = End
#Note: if the bin/fragment spans more than one 10kb region,
#then multiple 10 regions will be populated
sub add_bait {
    my ($csome, $start, $end, $hash_ref) = @_;
    
    my $ten_kb_start = ceil($start / 10_000);
    my $ten_kb_end = ceil($end/ 10_000);
    
    for (my $ten_kb_region = $ten_kb_start; $ten_kb_region <= $ten_kb_end; $ten_kb_region++){
        ${$hash_ref}{$csome}->{$ten_kb_region}->{$start} = $end;
    }
}



##########################################
#Subroutine: coord2bin
#Receives a chromosome name and a position and reference to the baits hash
#and returns the bait co-ordinates where this location is found (else returns 0)
#%lookup_hash{chromosome}->{10kb region}->{Start} = End
sub coord2bin{
    my ($csome, $pos, $hash_ref) = @_;
    my $ten_kb_region = ceil($pos / 10_000);

    foreach my $start ( keys %{ ${$hash_ref}{$csome}->{$ten_kb_region} }  ){
        my $end = ${ $hash_ref }{$csome}->{$ten_kb_region}->{$start};
        if ( ($start <= $pos) and ($end >= $pos) ){
            return ("$csome\t$start\t$end");
        }
    }
    return 0;    #Not captured
}


#Subroutine: calc_perc
#Receives a number and a total and returns the perentage value
#Optional argument: decimal places of the output
sub calc_perc {
    my ( $n, $total, $dp ) = @_;
    
  if(defined $dp){
      $dp = abs( int($dp) );
  }else{
      $dp = 2;
  }
  $dp = 1 / (10**$dp);   #'Nearest' function needs 1, 10, 100 etc. 
    
  return 'NA' unless(defined $n and defined $total);   #Avoid initialisation error
  return 'NA' if($total == 0);    #Avoid division by zero error
    
  my $pc = 100 * $n / $total;
  $pc = nearest($dp, $pc);
    
  return $pc;
}


#########################################
#Subroutine: deduplicate_array
#Takes and array and returns the array with duplicates removed
#(keeping 1 copy of each unique entry).
sub deduplicate_array {
    my @array = @_;
    my %uniques;

    foreach (@array) {
        $uniques{$_} = '';
    }
    my @uniques_array = keys %uniques;
    return @uniques_array;
}



1