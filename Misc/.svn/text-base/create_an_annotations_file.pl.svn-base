#!/usr/bin/perl

#############################################
#A Perl script to parse the file human_repeats_info_sequence.txt
#and remove selected repeat co-ordinvates
#
#The repeats file is of the format:
##bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd	repLeft	id
#0	1892	83	59	14	chr1	67108753	67109046	-181847376	+	L1P5	LINE	L15301	5607	-544	1
#1	2582	27	0	23	chr1	8388315	8388618	-240567804	-	AluY	SINE	Alu	-15	291
#1	4085	171	77	36	chr1	25165803	25166380	-223790042	+	L1MB5	LINE	L15567	6174	0	4
#
#Our annotations to process is of the format:
#repName	repClass	repFamily
#	DNA	
#	LINE	L1
#		L2
#		CR1
#Only one column is populated, the column to identify the Anaconda class
#
#
#The output file wil be tab-delimited, of the format:
#Chromosome    Start    End    Repeat_Name    Repeat_Class    Repeat_Family    Anaconda_Name   
#The Anaconda class will take its name from the input file (i.e. is EITHER the
#Name, Class or Family)
#

#############################################
use strict;
use warnings;

use Data::Dumper;

die "Please specify 2 files to process (1: list of repeats and 2: repeats of interest).\n" unless(scalar(@ARGV) == 2);

my ($repeats_file, $interest_file) = @ARGV;
my %feature_counter;    #%{feature} = count;


#Read Interest file
if ($interest_file =~ /\.gz$/ ) {
	open( INTEREST, "zcat $interest_file |" ) or die "Couldn't read file '$interest_file' : $!";
} else {
	open(INTEREST, '<', $interest_file) or die "Could not open '$interest_file' : $!";
}

print "Processing '$interest_file' to identify regions of interest\n";


my %interest_name;    # %{name} = ''
my %interest_class;    # %{class} = ''
my %interest_family;    # %{family} = ''

scalar <INTEREST>;    #Skip header
while(<INTEREST>){
	my $line = $_;
	chomp $line;

	my ($name, $class, $family) = split(/\t/, $line);
	$interest_name{$name} = '' if( (defined $name) and ($name ne '') );
	$interest_class{$class} = '' if( (defined $class) and ($class ne '') );
	$interest_family{$family} = '' if( (defined $family) and ($family ne '') );
}
close INTEREST or die "Could not close '$interest_file' : $!";


#print Dumper \%interest_name;
#print Dumper \%interest_class;
#print Dumper \%interest_family;


#Process repeats file that needs processing
print "Processing '$repeats_file' to identify interesting repeats\n";

if($repeats_file =~ /\.gz$/ ) {
	open( REPEATS, "zcat $repeats_file |" ) or die "Couldn't read file '$repeats_file' : $!";
} else {
	open( REPEATS, '<', $repeats_file) or die "Could not open '$repeats_file' : $!";
}

my $outfile = "$repeats_file.interest.txt.gz";
my $rejects_file = "$repeats_file.rejected.txt.gz";
open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
open(REJECTS, "| gzip -c - > $rejects_file") or die "Could not write to '$rejects_file' : $!";
print OUT "Chromosome\tStart\tEnd\tName\tClass\tFamily\tAnaconda_Name\n";


scalar <REPEATS>;    #Skip header
while(<REPEATS>){
	my $line = $_;
	chomp $line;

	my $csome = ( split(/\t/, $line) )[5];
	$csome =~ s/^chr//;
	my $start = ( split(/\t/, $line) )[6];
	my $end = ( split(/\t/, $line) )[7];
	my $name = ( split(/\t/, $line) )[10];
	my $class = ( split(/\t/, $line) )[11];
	my $family = ( split(/\t/, $line) )[12];

	if(exists $interest_name{$name}){
		print OUT "$csome\t$start\t$end\t$name\t$class\t$family\t$name\n";
	}elsif(exists $interest_class{$class}){
		print OUT "$csome\t$start\t$end\t$name\t$class\t$family\t$class\n";
	}elsif( exists $interest_family{$family} ){
		print OUT "$csome\t$start\t$end\t$name\t$class\t$family\t$family\n";
	}else{
		print REJECTS "$csome\t$start\t$end\t$name\t$class\t$family\t$family\n";
	}

}

close OUT or die "Could not close $outfile : $!";
close REJECTS or die "Could not close $rejects_file : $!";
print "Processing complete.\n";

exit (0);

