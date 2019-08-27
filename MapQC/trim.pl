#!/usr/bin/perl

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
