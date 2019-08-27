#!/usr/bin/perl

use strict;
use warnings;
use IO::File;                     # 5.004 or higher

use Data::Dumper;


unless(scalar @ARGV == 1){
        die "Please specify files to process.\n";
}

my $file = $ARGV[0];

if($file =~ /\.gz$/){
	open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
} else {
	open(IN, '<', $file) or die "Could not open '$file' : $!";
}

print "Processing $file\n";

my %tf_filehandle;	
while(<IN>){
	my $line = $_;
	chomp $line;
	my ($csome, $start, $end, $tf) = split(/\t/, $line);
	#($tf) = split(/_|\(/, $tf);
	
	if(exists $tf_filehandle{$tf}){    #Write to existing file   
		my $fh = $tf_filehandle{$tf}; 
		print $fh "$csome\t$start\t$end\t$tf\n";	
	}else{    #Write to new file    
		my $fh = IO::File->new();
		$tf_filehandle{$tf} = $fh;
		my $outfile = "$file.separate.$tf.bed.gz";
		open($fh, "| gzip -c - > $outfile") or die "Couldn't write to '$outfile' : $!";
		print $fh "$csome\t$start\t$end\t$tf\n";		
	}
}

print "Processing complete\n";

exit (0);