#!/usr/bin/perl

###########################################################
#Takes the fastq files to map, basenames of genome indices 
#and the splice site file and maps against the reference 
#genome using HISAT2
###########################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;

#Option variables
my %config = (
    genome => '',
    splice_sites => ''
);

my $config_result = GetOptions(
    "genome=s"   => \$config{genome},
    "splice=s"    => \$config{splice_sites},
);
die "Could not parse options.\n" unless ($config_result);

#Check input
unless(@ARGV){
	die "Please specify file(s) to process.\n";
}
die "Please adjust configuration.\n" unless ( check_files_exist(\@ARGV, 'EXISTS') );
#$config{genome} = '/bi/scratch/Genomes/Human/GRCh38/Homo_sapiens.GRCh38' unless( hasval($config{genome}) );
#$config{splice_sites} = '/bi/scratch/Genomes/Human/GRCh38/Homo_sapiens.GRCh38.83.splice_sites.txt' unless( hasval($config{splice_sites}) );

if($config{splice_sites} eq ''){
	die "Please specify a --splice file";
} else{
	unless(-e $config{splice_sites}){
		die "The --splice file '$config{splice_sites}' does not exist\n";
	}
}

	
foreach my $file (@ARGV){
	my $summary_filename = "$file.anaconda_mapper_summary_file.txt";
	open(SUMMARY, '>', $summary_filename) or die "Could not open '$summary_filename' : $!";
	
	#Perform mapping
	print "Mapping '$file' against $config{genome} with splice sites '$config{splice_sites}'\n";   
	my $outfile = $file;
    $outfile =~ s/\.fastq$//;
    $outfile .= '.sam';
	my $command = "hisat2 -x $config{genome} --known-splicesite-infile $config{splice_sites} -U $file -S $outfile 2>$summary_filename";
	print "Mapping with command: '$command'\n";
	!system($command) or die "Could not run command '$command'.\n";

    close SUMMARY or die "Could not close '$summary_filename' : $!";
}

print "Mapping complete\n";

exit (0);
