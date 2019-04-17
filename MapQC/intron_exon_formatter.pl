#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;

#######################################################################################
#Format the simulation formatted files that converts barcodes to numeric
#ids and edits the feature name with the type of feature (e.g. NASCENT) and then creates final
#simulation output files
#
#The script opens the input files and then relates all sequence barcodes to a numerical index
#The files are then re-written and features are suffixed with an identifier e.g. NASCENT

#######################################################################################


########################################################
#Get user-supplied parameters
#Option variables
my %config = (
    outfile    => undef,
);

my $config_result = GetOptions(    #Stores parameters
    "outfile=s"   => \$config{outfile},
);

die "Command line options need to be in the correct format.\n" unless ($config_result);


#Check input ok
unless (@ARGV){
    die "Please specify files to process.\n";
}
my @files = deduplicate_array(@ARGV);


die "Specify an output file.\n" unless defined($config{outfile});


#Process files
my %barcodes;    # %{SEQUENCE} = index
my $id = 0;


open(OUT, "| gzip -c - > $config{outfile}") or die "Could not write to '$config{outfile}' : $!";
my $header;

print "Processing file:\n";

foreach my $file (@files){
	print "\t$file\n";
	my $suffix;
	if($file =~ /EXON_JUNCTION/){
		$suffix = 'EXON_JUNCTION';
	} elsif($file =~ /NASCENT/){
		$suffix = 'NASCENT';
	} elsif($file =~'PUTATIVE_EXON'){
		$suffix = 'PUTATIVE_EXON';
	} elsif($file =~ 'OTHER'){
		$suffix = 'OTHER';
	} else {
		die "Could not fild valid suffix in '$file'.\n";
	}

	

#EXON_JUNCTION => 0, NASCENT => 0, PUTATIVE_EXON => 0, OTHER 

	my $fh_in = cleverOpen($file); 

	if(!defined $header){    #Print header once only
		$header = scalar <$fh_in>;
		print OUT $header
	} else {
		scalar <$fh_in>;
	}

	while(<$fh_in>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		my $barcode = $line_elements[1];
		my $id_to_use;

		if(exists $barcodes{$barcode}){
			$id_to_use = $barcodes{$barcode};
		} else{
			$id++;	
			$barcodes{$barcode} = $id;
			$id_to_use = $id;
		}
	
		$line_elements[1]  = $id_to_use;
		$line_elements[7] = "$line_elements[7].$suffix";
		print OUT join("\t", @line_elements) . "\n"; 
	}

	close $fh_in or die "Could not close '$file' filehandle : $!";
}

close OUT or die "Could not close filehandle on '$config{outfile}' : $!";

print "Created formated exon/intron simulation files.\n";

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

