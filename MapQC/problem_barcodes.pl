#!/usr/bin/perl

##############################################################
#Takes a FASTQ file in which the first 26bps constitute
#a barcode. The (trimmed i.e. 20bps) barcode sequence is then extracted
#and sequences with an extreme A:G:C:T proportion are removed.
#
#Also maps the barcode against possible sources of contamination
#e.g. adapter sequences.
#
#
##############################################################

use strict;
use warnings;
use String::Approx 'amatch';
use Math::Round;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;

unless(@ARGV){
	die "Please specify a file to process.\n";
}

my %summary_counter;
my @summary_categories = qw( Reads_Processed Fail_multiple_bases Fail_problem_sequence
						Pass Percent_Pass);


print "Checking for problem barcode sequences in...\n";
foreach my $file (@ARGV){

	print "\t$file\n";
	 
	if ( $file =~ /\.gz$/ ) {
		open( IN, "zcat $file |" ) or die "Couldn't read file '$file' : $!";
	} else {
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}
	
	my $pass_file = $file;
	$pass_file =~ s/\.present\.fastq\.gz$//;
	$pass_file .= '.barcode_pass.fastq.gz';
	
	my $fail_file = $file;
	$fail_file =~ s/\.present\.fastq\.gz$//;
	$fail_file .= '.barcode_fail.fastq.gz';
		
	
	my $summary_file = "$file.problem_barcodes_summary_file.txt";
	open( PASS, "| gzip -c - > $pass_file" ) or die "Couldn't write to file '$pass_file': $!";
	open( FAIL, "| gzip -c - > $fail_file" ) or die "Couldn't write to file '$fail_file': $!";
	open( SUMMARY, '>', $summary_file) or die "Couldn't write to file '$summary_file' : $!";
	
	print SUMMARY "File";
	foreach my $category (@summary_categories){   #Reset summary counter and print headers
		$summary_counter{$category} = 0;
		print SUMMARY "\t$category";
	}
	print SUMMARY "\n";

	while(<IN>){
		my $read = $_;
		$read .= scalar <IN>;
		$read .= scalar <IN>;
		$read .= scalar <IN>;
			
		my $barcode = (split(/\n/, $read))[1];
		$barcode = substr($barcode, 0, 26);
		my $barcode_trimmed = substr($barcode, 3, 20);    #Trim the barcode prior to mapping

		$summary_counter{Reads_Processed}++;
	
		#Does one base occur more than expected?
		if(extreme_base($barcode_trimmed)){
			$summary_counter{Fail_multiple_bases}++;
			print FAIL $read;
			next;		
		}
	
		#Possible contaminants?
		if(contamination($barcode_trimmed)){
			$summary_counter{Fail_problem_sequence}++;
			print FAIL $read;
			next;
		}
		
		print PASS $read;
		$summary_counter{Pass}++;
		
	}
	close PASS or die "Could not close filehandle on '$pass_file' : $!";
	close FAIL or die "Could not close filehandle on '$fail_file' : $!";

	#Print summary results
	$summary_counter{Percent_Pass} = calc_perc($summary_counter{Pass}, $summary_counter{Reads_Processed}, 1);
	print SUMMARY "$file";
	foreach my $category (@summary_categories){
		print SUMMARY "\t$summary_counter{$category}";
	}
	print SUMMARY "\n";
	#print Dumper \%summary_counter;

	close SUMMARY or die "Could not close filehandle on '$summary_file' : $!";


}	

print "Checking complete.\n";

exit (0);




##########################################################
#Subroutines
##########################################################

#Subroutine: contamination
#Takes a sequence and checks whether it is a likely contaminant
#Returns 1 if it is a possible contaminant, 0 if not
sub contamination{
	my $barcode = $_[0];
	$barcode = substr($barcode, 2, 16);    #Trim a little for the lookup
	
	my @contaminants = qw( 
		CCATCTCATCCCTGCGTGTC
		GATCGTCGGACTGTAGAACTCCCTATAGTGAGTCGTATTACAAGGCACACAGGGGATAGG
		CCACTACGCTCGCTATCCTATCCCCTGTGTGCCTTG 
		CCATCTCATCCCTGCGTGTC 
		AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGA
		CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
	);
	
	foreach my $contaminant (@contaminants){
		return 1 if($contaminant =~ /$barcode/);
		my $reverse_complement_contaminant = revcom($contaminant);
		return 1 if($reverse_complement_contaminant=~ /$barcode/);
	}
	return 0;
}



#Subroutine: revcom
#Reverse complements a sequence
sub revcom{
		my $seq = $_[0];
		$seq = reverse($seq);
		$seq =~ tr/AGCTagct/TCGAtcga/;
		return $seq;
}



#Subroutine: extreme_base
#Takes a polynucleotide sequence and determines whether one base occurs
#more than expected by chance (predetermined threshold). If it does the
#subroutine returns 1, if not it returns 0.
sub extreme_base {
	my $seq = $_[0];
	my $threshold = 13;
	$seq = uc $seq;	
	my $max_count = 0;
	my $count;
	
	$count = ($seq =~ tr/A//);
	$max_count = $count if($count > $max_count);
	$count = ($seq =~ tr/C//);
	$max_count = $count if($count > $max_count);
	$count = ($seq =~ tr/G//);
	$max_count = $count if($count > $max_count);
	$count = ($seq =~ tr/T//);
	$max_count = $count if($count > $max_count);
	
	if($max_count > $threshold){
		return 1;
	}else{
		return 0;
	}
}
	