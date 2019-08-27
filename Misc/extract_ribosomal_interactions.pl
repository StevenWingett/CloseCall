#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

#######################################################################################
#Process an Anaconda format file to identify groups containing ribosomal containing features
#These groups are printed to an output file, minus the ribosomal features
#######################################################################################


my @desired_features = ('6.28S', '2.18S');


#Check input ok
unless (@ARGV){
    die "Please specify files to process.\n";
}
my @files = deduplicate_array(@ARGV);



#Process files
print "Processing:\n";
foreach my $file (@files){

	my %contains_rib;    # %{barcode_id} = '';   If contains ribosomal feature


	print "\t$file\n";
	
	#Pass 1
	print "\t\tIdentifing barcode groups\n";
	my $fh_in = cleverOpen($file); 
	scalar <$fh_in>;   #Ignore header

	while(<$fh_in>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		my $bacodeID = $line_elements[1];
		my $feature_name = $line_elements[-1];
		
		foreach my $desired_feature (@desired_features){
			if($desired_feature eq $feature_name){
				$contains_rib{$bacodeID} = '';
				last;
			}
		}
	}
	close $fh_in or die "Could not close '$file' filehandle : $!";
	
	

	#print Dumper \%contains_rib;
	
	
	#Pass 2
	print "\t\tIdentifing barcode groups\n";
	$fh_in = cleverOpen($file); 
	my $outfile = "$file.ribosomal.edited.txt";
	open(OUT, '>', $outfile) or die "Could not open filehandle on '$outfile' : $!";
	print OUT scalar <$fh_in>;   #Ignore header
	while(<$fh_in>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		my $bacodeID = $line_elements[1];
		next unless(exists $contains_rib{$bacodeID});
		
		my $feature_name = $line_elements[-1];
		
			
		my $print_flag = 1;	
		foreach my $desired_feature (@desired_features){
			if($desired_feature eq $feature_name){
				$print_flag = 0;

				
				
			}
		}
		
		print OUT "$line\n" if $print_flag;
	
	}
	
	
	close OUT or die "Could not close filehandle on '$outfile' : $!";
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
  
	if ($file =~ /\.gz$/){
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

