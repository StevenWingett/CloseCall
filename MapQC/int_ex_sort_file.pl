#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;

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


#######################################################################################
#Perl script to sort the simulation formatted file so it is in the correct order for the simulation
#Will report if duplicates are identified, since these should not be present
#######################################################################################


########################################################
#Get user-supplied parameters
#Option variables
my $outfile;

my $config_result = GetOptions(    #Stores parameters
    "outfile=s"   => \$outfile,
);

die "Command line options need to be in the correct format.\n" unless ($config_result);

die "Specify an output file.\n" unless defined($outfile);


#Check input ok
unless (@ARGV){
	die "Please specify files to process.\n";
}
my @files = deduplicate_array(@ARGV);





#Process files
print "Processing:\n";
foreach my $file (@files){

	print "\t$file\n";
	my $fh_in = cleverOpen($file); 

	open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
	print OUT scalar <$fh_in>;   #Header

	my %data;   # %{barcode_ID} = @data
	
	while(<$fh_in>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		shift @line_elements;   #Remove read ID
		my $barcode_ID = shift @line_elements;
		my $info = join("\t", @line_elements);
		
		push (@{ $data{$barcode_ID} }, $info);

	}
	close $fh_in or die "Could not close '$file' filehandle : $!";
	
	#Sort data	
	my $readID = 1;
	
	my %dedup_check;   #%{barcode\tfeature};
	my $dups = 0;
	
	foreach my $barcode_ID (sort {$a <=> $b} keys %data){
		foreach my $info ( @{ $data{$barcode_ID} } ){
			my $feature = (split(/\t/, $info))[5];
			my $lookup = "$barcode_ID\t$feature";
			#print "$lookup\n";
			if(exists $dedup_check{$lookup}){
				#print "$readID\t$barcode_ID\t$info\n";
				$dups++;
			} else{
				print OUT "$readID\t$barcode_ID\t$info\n";
				$readID++;
				$dedup_check{$lookup} = '';
			}
		}
	}
	
	close OUT or die "Could not close filehandle on '$outfile' : $!";
	
	#print Dumper \%data;

print "Warning!  Duplicates identified\n" if($dups);
print "Number of duplicates: $dups\n";

}

print "Processing complete.\n";

exit (0);


#####################################################################
#Subroutines
#####################################################################

#######################
##Subroutine "cleverOpen":
##Opens a file with a filhandle suitable for the file extension
sub cleverOpen{
  my $file  = shift;
  my $fh;
  
	if( $file =~ /\.bam$/){
		open( $fh, "samtools view -h $file |" ) or die "Couldn't read '$file' : $!";  
	}elsif ($file =~ /\.gz$/){
		open ($fh,"zcat $file |") or die "Couldn't read $file : $!";
	} else {
		open ($fh, $file) or die "Could not read $file: $!";
	}
  return $fh;
}



#Sub: deduplicate_array
#Takes and array and returns the array with duplicates removed
#(keeping 1 copy of each unique entry).
sub deduplicate_array{
	my @array = @_;
	my %uniques;

	foreach (@array){
		$uniques{$_} = '';	
	}
	my @uniques_array = keys %uniques;
	return @uniques_array;
}

