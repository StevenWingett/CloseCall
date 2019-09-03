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


use strict;
use warnings;

######################################################
#Script to remove the first 67 bp from the start of a 
#read and keep the next 50 bps and remove the rest 
######################################################

#Pass file names as command-line arguments
my @files = @ARGV;

foreach my $file (@files){
  
  chomp $file;
  
  warn "Processing $file\n";
  
  if ($file =~ /\.gz$/){
    open (IN,"zcat $file |") or die "Couldn't read $file : $!";  
  }else{
    open (IN, $file) or die "Could not open $file\n";
  }

  my $outfile = $file;
  $outfile =~ s/\.barcode_pass\.fastq\.gz\.barcoded\.fastq$//;
  $outfile .= '.trimmed.fastq';
  open(OUT, '>', $outfile) or die "Couldn't write to '$outfile' : $!";

  while(<IN>){
    if(/^@/){
      my $line1 = $_;
      my $line2 = scalar <IN>;
      my $line3 = scalar <IN>;
      my $line4 = scalar <IN>;

      
      $line2 = substr($line2, 67, 50);
      $line4 = substr($line4, 67, 50);

      print OUT $line1 . $line2. "\n" . $line3 . $line4."\n";

    }else{
      die "@ not found at start of read!\n";
    }
  }

  close IN;
  close OUT or die "Could not close filehandle out '$outfile' : $!";

}
 
print "Processing compete\n";

exit (0);
