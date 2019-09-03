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


######################################################
#Takes basenames of genome indices and maps against
#these. When mapping against the first genome, Anaconda 
#does not allow multi-mapping reads, but does against the
#remaining genomes. The output files are then processed to
#allow only reads which map to one genome.
#Takes a list of files to process and genome basenames:
#Genome1 Genome2 Genome3
#The first genome in the list is the one in which multi-mapping
#is not allowed
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib $Bin;
use File::Basename;

use Data::Dumper;

#Option variables
my %config = (
    genome => '',
    repeats   => '',
);

my $config_result = GetOptions(
    "genome=s"   => \$config{genome},
    "repeats=s"    => \$config{repeats},
);
die "Could not parse options.\n" unless ($config_result);


#Check input
unless(@ARGV){
	die "Please specify file(s) to process.\n";
}
die "Please adjust configuration.\n" unless ( check_files_exist(\@ARGV, 'EXISTS') );
$config{genome} = '/bi/scratch/Genomes/Human/GRCh38/Homo_sapiens.GRCh38' unless( hasval($config{genome}) );
my @repeats;
@repeats = split(/,/, $config{repeats}) if( hasval($config{repeats}) );




	
foreach my $file (@ARGV){

	my $summary_filename = "$file.anaconda_mapper_summary_file.txt";
	open(SUMMARY, '>', $summary_filename) or die "Could not open '$summary_filename' : $!";
	print SUMMARY "File\tReads_Mapping_1_Genome\tReads_Mapping_Multiple_Genomes\n";

	my @outfiles;
	#Perform mapping
	print "Mapping '$file' against:\n";
	print "\t$config{genome}\n";    #Multi-mapping NOT allowed
	my $outfile = basename($file);
	$outfile = "$outfile." . basename($config{genome}) . ".sam";
	push (@outfiles, $outfile);
	my $command = "bowtie -m 1 -p 4 --best --sam --chunkmbs 512 $config{genome} $file $outfile >/dev/null";
	!system($command) or die "Could not run command '$command'.\n";

	foreach my $repeat_genome (@repeats){
		print "\t$repeat_genome\n";
		my $outfile = basename($file);
		$outfile = "$outfile." . basename($repeat_genome) . ".sam";
		push (@outfiles, $outfile);
		my $command = "bowtie -k 1 -p 4 --best --sam --chunkmbs 512 $repeat_genome $file $outfile >/dev/null";
		!system($command) or die "Could not run command '$command'.\n";
	}


	#Identify reads in output files
    my %read_counter;    #%{id} = count
    foreach my $outfile (@outfiles){ 

        if ( $outfile =~ /\.gz$/ ) {
            open( IN1, "zcat $outfile|" ) or die "Couldn't read $outfile : $!";
        } elsif ( $outfile =~ /\.bam$/ ) {
             open( IN1, "$config{samtools} view -h $file |" ) or die "Couldn't read $file: $!";
        } else {
            open( IN1, $outfile ) or die "Could not read $outfile : $!";
        }

        while(<IN1>){
            my $line = $_;
            chomp $line;
            next if(substr($line, 0, 1) ) eq '@';
            my ($id, $samflag) = split(/\t/, $line);

            next if ($samflag & 0x4);
            $read_counter{$id}++;
        }

        close IN1 or die "Could not close filehandle on '$outfile'\n";
    }

	#Write reads that map to one genome to an output file

    my $combined_outfile = "$file.mapped.sam";
    my %results = ('Reads_Mapping_1_Genome' => 0, 'Reads_Mapping_Mutiple_Genomes' => 0);
	
	
	#Print headers to final outputfile
	open(OUT, '>', $combined_outfile) or die "Could not write to '$combined_outfile'.\n";
    foreach my $outfile (@outfiles){ 
		if ( $outfile =~ /\.gz$/ ) {
            open( IN2, "zcat $outfile|" ) or die "Couldn't read $outfile : $!";
        } elsif ( $outfile =~ /\.bam$/ ) {
            open( IN2, "$config{samtools} view -h $file |" ) or die "Couldn't read $file: $!";
        } else {
            open( IN2, $outfile ) or die "Could not read $outfile : $!";
        }
	
	    while(<IN2>){
			my $line = $_;
			chomp $line;
			
			if(substr($line, 0, 1) eq '@'){    #Header
				print OUT "$line\n" ;
				last;
			}
		}
	}
	close IN2 or die "Could not close filehandle on '$outfile'\n";
	
	
	#Print out reads mapping to only one genome to outputfile
    open(OUT, '>', $combined_outfile) or die "Could not write to '$combined_outfile'.\n";
    foreach my $outfile (@outfiles){ 

        if ( $outfile =~ /\.gz$/ ) {
            open( IN3, "zcat $outfile|" ) or die "Couldn't read $outfile : $!";
        } elsif ( $outfile =~ /\.bam$/ ) {
             open( IN3, "$config{samtools} view -h $file |" ) or die "Couldn't read $file: $!";
        } else {
            open( IN3, $outfile ) or die "Could not read $outfile : $!";
        }

        while(<IN3>){
            my $line = $_;
            chomp $line;
			
			if(substr($line, 0, 1) eq '@'){    #Header
				print OUT "$line\n" ;
				next;
			}
			
            my ($id, $samflag) = split(/\t/, $line);

            next if ($samflag & 0x4);
            if(exists $read_counter{$id}){
                if($read_counter{$id} == 1){
                    print OUT "$line\n";
					$results{Reads_Mapping_1_Genome}++;
                }elsif($read_counter{$id} > 1)  {
					$results{Reads_Mapping_Mutiple_Genomes}++;
                }
            }else{
				die "'$read_counter{$id}' not found in lookup hash!\n";
			}			
        }
        close IN3 or die "Could not close filehandle on '$outfile'\n";
    }

    close OUT or warn "Could not close filehandle on '$combined_outfile'\n";
	print SUMMARY "$file\t$results{Reads_Mapping_1_Genome}\t$results{Reads_Mapping_Mutiple_Genomes}\n";
	
	close SUMMARY or die "Could not close '$summary_filename' : $!";
}



print "Mapping complete\n";

exit (0);














############################################################
#Subroutines
############################################################


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


