#!/usr/bin/perl

#############################################
#Append 1 billion to fist column in file
#############################################

use strict;
use warnings;

use Data::Dumper;

unless(@ARGV){
	die "Please specify a file to process.\n";
}

foreach my $file(@ARGV){
	my $outfile = "appended.$file";
	open(IN, '<', $file) or die "Could not open '$file' : $!";
	open(OUT, '>', $outfile) or die "Could not write to '$outfile' : $!";

	while(<IN>){
		my $line = $_;
		chomp $line;
		my @data = split(/\t/, $line);
		$data[0] = $data[0] + 1_000_000_000;
		print OUT join("\t", @data);
		print OUT "\n";
	}

	
	close IN or die "Could not close $file : $!";
	close OUT or die "Could not close $outfile : $!";


}


print "Processing complete\n";

exit (0);
