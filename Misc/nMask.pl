#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;

###################################################################################
#Perl script to N-mask FASTA files
##################################################################################

#Get user-supplied parameters
#Option variables
my %config = (
    help   => '',
    regions => ''
);

my $config_result = GetOptions(
    "help"    => \$config{help},
    "regions=s" => \$config{regions}
);

die "Command line options need to be in the correct format.\n" unless ($config_result);


if ( $config{help} ) {
    print while (<DATA>);
    exit;
}

unless(@ARGV){
	die "Please specify file(s) to process.\n";
}

die "Specify a --regions file.\n" if ($config{regions} eq ''); 

#Extract data from the FASTA files
my %chromosomes_processed;    #%{Chromosome_Name} = sequence
foreach my $file (@ARGV) {

	print "Reading '$file'\n";

    if ( $file =~ /.*\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Cannot open filename: $!";
    } else {
        open( IN, $file ) or die "Can't read filename: $!";
    }
	
    my $chromosome = '';
    while (<IN>) {
        my $line = $_;
        chomp($line);

        if ( $line =~ /^\>/ ) {    #Process FASTA header to determine chromosome
        	unless($line =~ /^\>(\S+).*$/){
        		die "File '$file' is not in FASTA format at line:\n$line\nA chromosome name is reqired immediately after '>' (blank spaces are not allowed)\nFor more details on FASTA format go to:\nhttp://en.wikipedia.org/wiki/FASTA_format\n";
        	}
			
			$chromosome = $1;
			
			if ( exists $chromosomes_processed{$chromosome} ) {
				die "The sequence file(s) to be digested contain multiple instances of chromosome $chromosome. Digestion terminated.\n";
            } else {
				print "\tFound chromsome: '$chromosome'\n";
                $chromosomes_processed{$chromosome} = '';
            }
			
        } else {
			die "Chromsome not defined in FASTA header line in file '$file'\n" if($chromosome eq '');
            $chromosomes_processed{$chromosome} .= $line;
        }
    }
    close IN or die "Could not close '$file' : $!";
}

#Determine regions to N-mask
print "Reading regions file $config{regions}\n";

if ( $config{regions} =~ /.*\.gz$/ ) {
    open( REGIONS, "zcat $config{regions} |" ) or die "Cannot open '$config{regions}' : $!";
} else {
    open( REGIONS, $config{regions}) or die "Can't read '$config{regions}' : $!";
}

scalar <REGIONS>;    #Ignore header
while(<REGIONS>){
	my $line = $_;
	chomp $line;
	next if($line =~ /^\s*$/);
	my ($csome, $start, $end) = split(/\t/, $line);
	my $length = $end - $start + 1;
	my $poly_n = 'N' x $length;
	substr($chromosomes_processed{$csome}, ($start - 1), $length, $poly_n);    #Note substr usage (assignment operator not needed)
}

close REGIONS or die "Could not close '$config{regions}' : $!";


#Write out resuls
my $outfile = 'nMasked_genome.fasta';
print "Writing results to $outfile\n";
open(MASKED, '>', $outfile) or die "Could not write to '$outfile' : $!";
#my $fasta_length = 6;   #Base pairs printed on single line in output FASTA file
foreach my $csome (sort keys %chromosomes_processed){
	print "\tWriting masked chromsome '$csome'\n";
	print MASKED ">$csome\n";    #Header
	my @split_data = unpack("(a100)*", $chromosomes_processed{$csome});
	foreach my $line(@split_data){
		print MASKED "$line\n";
	}
}

close MASKED or die "Could not close '$outfile' : $!";
print "Processing complete.\n";

exit (0);


__DATA__

FUNCTION

Perl Script to N-mask specified regions of the genome

SYNOPSIS

nMask.pl --regions [FASTA files]

COMMAND LINE OPTIONS

--regions        Tab-delimied file of positions to N-mask
                 Chromosome Start End (1-based co-ordinates)
--help           Print help message and exit

Steven Wingett, Babraham Institute, Cambridge, UK (steven.wingett@babraham.ac.uk)
