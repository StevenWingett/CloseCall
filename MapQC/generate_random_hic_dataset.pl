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


######################################################################
#Takes genome FASTA files (to determine chromsome lenghts and creates
#a virtual HiC dataset
######################################################################

use strict;
use warnings;
use Getopt::Long;
use Math::Random;
use POSIX;    #Needed only for custom script

use Data::Dumper;

#Get user-supplied parameters
#Option variables
my %config = (
    ditags   => '',
    help	=> ''
);

my $config_result = GetOptions(
    "ditags=i"    => \$config{ditags},
    "help"	=> \$config{help}
);

die "Command line options need to be in the correct format\n" unless ($config_result);

if ( $config{help} ) {
    print while (<DATA>);
    exit(0);
}

die "Please specify genome FASTA files fo process.\n" unless(@ARGV);

if($config{ditags} eq ''){
		print "Number of --ditags defauting to 100 million\n";
		$config{ditags} = 100_000_000;
}

die "Number of ditags needs to be at least 1.\n" if($config{ditags} < 1);

#Read FASTA files to determnine chromosome lenghts
my %csome_length = get_csome_lengths(@ARGV);    #%{chromosome} = length (bps)
#%csome_length = modify_for_long_csomes(%csome_length);







#Create the random ditags file
my $outfile = "$config{ditags}.random_ditags.txt";
generate_random_reads(\%csome_length, ($config{ditags} * 2), $outfile );    #ditags = reads x 2

print "Processing complete.\n";

exit (0);






##############################################################################################
#Subroutines
##############################################################################################

#Needed?
# Subroutine modify_for_long_csomes
# Takes  %csome_length and subdivides chromsome names into 2_000_000_000 chunks, since
# the random number gernerator has a limit of 2147483561
# e.g. {X} = 2_000_000_001 
# Becomes {X\t1} = 2_000_000_000 and {X\t2} = 1
# sub modify_for_long_csomes{
	# my %orig_hash = @_;
	# my %new_hash;
	# my $chunk_size = 2_000_000_000;
	
	# foreach my $csome (keys %orig_hash){
		# my $length = $orig_hash{$csome};
		
		# my $position = 0;
		# my $i = 0;
		
		# do{		

			# $i++;
			# $position = $i * $chunk_size;
		
			# if( ($chunk_size * $i) > $length ){
				# my $final_length = $length - ($chunk_size * ($i - 1));    #This length is smaller than the chunk size
				# $new_hash{"$csome\t$i"} = $final_length;
			# }else{
				# $new_hash{"$csome\t$i"} = $chunk_size;
			# }


		# }until($position >= $length)
	# }
	# return %new_hash;
# }


#Subroutine  generate_random_reads:
sub generate_random_reads{

 	my ( $csome_length_ref, $reads, $outfile) = @_;
 	my $reads_generated = 0;

 	my $start = 1;
 	my $end;
 	my %csome_running_size;    #Adds chromsome lengths sequentially e.g. X => 1..100, Y => 101..200
 	foreach my $csome (sort keys %{ $csome_length_ref} ){
		$end = $start + ${ $csome_length_ref }{$csome} - 1;
		$csome_running_size{"$start\t$end"} = $csome;
 		$start = $end + 1;
 	}
 	my $total_genome_size = $end;

 	print "Creating $reads random positions\n";
	
	my @random_positions = random_uniform_integer($reads, 1, int($total_genome_size / 2)  );    #Divide by 2 so not exceed limit for very large genomes
	foreach my $position (@random_positions){
		$position = $position * 2;
	}
	
	#Convert positions into a chromosome/position
	print "Converting random reads to genomic positions\n";
	open( OUT, '>', $outfile) or die "Could not write to '$outfile' : $!";
	foreach my $random_position (@random_positions){
		$reads_generated++;
		unless($reads_generated % 10_000_000){
			print "$reads_generated reads converted to genomic positions\n";
		}

		foreach my $position (keys %csome_running_size){	
			my ($start, $end) = split(/\t/, $position);
			if( ($random_position >= $start) and ($random_position <= $end) ){
				my $csome = $csome_running_size{$position};
				my $csome_position = $random_position - $start + 1;
#######			print OUT "$csome\t$csome_position\n";    #This should be standard
				
				my $group = ceil($reads_generated / 2);
				print OUT "$group\t$csome\t$csome_position\t+\n";    #Customised for later script 
			}
		}
	}
	close OUT or die "Could not close filehandle on '$outfile' : $!";
}



#Subroutine get_csome_lengths:
#Takse a list of FASTA files and returns as hash:
#%{chromosome_name} = chromosome_length
sub get_csome_lengths{
	my @filenames = @_;	
	my %chromosomes_processed;
	
	foreach my $file (@filenames) {

		if ( $file =~ /.*\.gz$/ ) {
			open( IN, "zcat $file |" ) or die "Cannot open file: $!";
		} else {
			open( IN, $file ) or die "Can't read file: $!";
		}
		
		warn "Reading '$file'\n";

		my $sequence = '';
		my $chromosomes_in_file = 0;
		my $line_count = 0;
		my $csome;
		
		while (<IN>) {
		
			my $line = $_;
			chomp($line);

			$line_count++;

			if($line_count == 1){    #Check first line is a FASTA header line
				unless($line =~ /^\>(\S+).*$/){
					die "File '$file' is not in FASTA format at line:\n$line\nFirst line of a FASTA file should give the chromosome name\nFor more details on FASTA format go to:\nhttp://en.wikipedia.org/wiki/FASTA_format\n";
				}
			}

			if ( $line =~ /^\>/ ) {    #Process FASTA header to determine chromosome
				unless($line =~ /^\>(\S+).*$/){
					die "File '$file' is not in FASTA format at line:\n$line\nA chromosome name is reqired immediately after '>' (blank spaces are not allowed)\nFor more details on FASTA format go to:\nhttp://en.wikipedia.org/wiki/FASTA_format\n";
				}

				if ( $chromosomes_in_file > 0 ) {    #More than 1 c'some in file - process previous c'some now
					$chromosomes_processed{$csome} = $sequence;
					$sequence = '';
				}

				#Current c'some identity
				$csome = $1;
				$chromosomes_in_file++;
				if ( exists $chromosomes_processed{$csome} ) {
					die "The sequence file(s) to be digested contain multiple instances of chromosome '$csome'. Digestion terminated.\n";
				} else {
					$chromosomes_processed{$csome} = '';
				}

			} else {
				$sequence .= $line;
			}
		}
		close IN;
		$chromosomes_processed{$csome} = $sequence;    #The last (or only) sequence in the file
	}
	
	#Determine lengths
	my %csome_length;
	foreach my $csome (keys %chromosomes_processed){
		$csome_length{$csome} = length ( $chromosomes_processed{$csome} );
	}
	return %csome_length;
}
		



__DATA__

SYNOPSIS

Creates random Hi-C ditags

generate_random_hic_dataset [OPTIONS]... [FASTA FILES]

FUNCTION

Uses FASTA files to determine the length of each chromosome and then generates a
random Hi-C dataset. Duplicate di-tags may be generated.

COMMAND LINE OPTIONS

--ditags       Number of ditags to create in Hi-C dataset
--help         Print help message and exit

Steven Wingett, Babraham Institute, Cambridge, UK (steven.wingett@babraham.ac.uk)

