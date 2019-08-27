#!/usr/bin/perl

##############################################################################
#Perl script takes a read and extracts the barcode (in-line, 26bps) from
#the read (single-end) and places it in the header of the read. 
#This script uses another file which in which each barcoded has been assigned
#a barcode group - this is what is put in the header.
#Format:
#Barcode \t Group
#The name of this barcode->group conversion file is deduced from the sequence
#input files
##############################################################################

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;

#Check input ok
die "Please specify files to process to process.\n" unless (@ARGV);
die "Please adjust configuration.\n" unless( check_files_exist(\@ARGV, 'EXISTS') );


#Append the barcode from ReadF to ReadR
foreach my $file (@ARGV){

	#Create summary file
	my @summary_categories = qw(File Reads_Processed Reads_assigned_barcode Reads_not_assigned_barcode);
	open(SUMMARY, '>', "$file.assigned_barcodes_summary.txt") or die "Could not write to '$file.assigned_barcodes_summary.txt' : $!";
	my $summary_header = '';
	foreach my $category ( @summary_categories  ){
		$summary_header .= "$category\t";
	}
	$summary_header =~ s/\t$/\n/;    #Replace last tab with new line
	print SUMMARY $summary_header;

	my %summary_counter;
	foreach my $category (@summary_categories){
		$summary_counter{$category} = 0;
	}
		
	print "Appending barcodes in $file to sequences in $file\n";
	
	if ( $file=~ /.*\.gz$/ ) {
        open( FASTQ_IN, "zcat $file |" ) or die "Could not open $file : $!";
    } else {
        open( FASTQ_IN, $file ) or die "Could not open $file : $!";
    }
	
	$summary_counter{File} = $file;
	my $outFilename = $file . '.barcoded.fastq';

    #Open Barcode -> group conversion file and determine conversions
    my $conversion_file = $file;
    $conversion_file =~ s/\.barcode_pass\.fastq\.gz$/\.barcode_list\.fasta/;
   # $conversion_file =~ s/\.barcoded.fastq$//;
    $conversion_file .= '.sam.barcode_groups.txt';
    print "Reading conversion file '$conversion_file'\n"; 
    open( CONVERSION, '<', $conversion_file ) or die "Could not open '$conversion_file' : $!";
    my %barcode_group;
    while(<CONVERSION>){
        chomp;
        my ($barcode, $group) = split(/\t/);
        if(exists $barcode_group{$barcode}){
            die "Error: '$barcode_group{$barcode} does not match '$group'.\n" if($barcode_group{$barcode} ne $group);    #Check input ok
        }
        $barcode_group{$barcode} = $group;
    }
    close CONVERSION or die "Could not close filehandle on '$conversion_file' : $!";
	
	open( OUT, '>', $outFilename ) or die "Couldn't read file '$outFilename : $!";
	
	while(<FASTQ_IN>){

		my $line1 = $_;
		my $line2 = scalar <FASTQ_IN>;
		my $line3 = scalar <FASTQ_IN>;
		my $line4 = scalar <FASTQ_IN>;
		
		$summary_counter{Reads_Processed}++;
				
		my $barcode = substr($line2, 0, 26);    #Line2
		

        if(exists $barcode_group{$barcode}){
			$summary_counter{Reads_assigned_barcode}++;
            $barcode = $barcode_group{$barcode};    #Convert to barcode group (only those assigned a group)
        }else{
			$summary_counter{Reads_not_assigned_barcode}++;
            next;    #Move to next read
        }

		$line1 = substr($line1, 1);    #Remove starting @
		$line1 = '@' . $barcode . '_' . $line1;   #Add barcode at start, and then prefix with '@' syymbol
			
		print OUT $line1.$line2.$line3.$line4 if $barcode =~ /^[ATGC]+$/;    #Ignore barcodes containing 'N'
			
	}
	
	close FASTQ_IN or die "Could not close filehandle on '$file' : $!";
	close OUT or die "Could not close filehandle on '$outFilename' : $!";
	
	#Print results to summary file
	my $summary_line = '';
	foreach my $category (@summary_categories ){
		$summary_line .= "$summary_counter{$category}\t";
	}
	$summary_line =~ s/\t$/\n/;    #Replace last tab with new line
	print SUMMARY $summary_line;	
	close SUMMARY or die "Could not close filhandle on '$file.assigned_barcodes_summary.txt' : $!";
}

print "Processing complete.\n";

exit (0);



#######################################################################################
#Subroutines                                                                          #
#######################################################################################


#############################################################
#check_files_exist:
#Takes a reference to an array containing paths to filenames and verifies they exist
#Warns of files that do no exit. Returns 1 if all files exist but 0 if this is not
#the case.
#
#Also, takes a second argument:
#$_[1] should be 'EXISTS' or 'NOT_EXISTS'
#If 'NOT_EXIST' warns if file already exists.  Returns '1' if none of the
#files exists and '0' if one or multiple files already exist
# sub check_files_exist {
#     my $files      = $_[0];    #Reference to array
#     my $check_for  = $_[1];
#     my $all_exist  = 1;
#     my $not_exists = 1;

#     if ( $check_for eq 'EXISTS' ) {
#         foreach my $file (@$files) {
#             unless ( -e $file ) {
#                 warn "$file does not exist\n";
#                 $all_exist = 0;
#             }
#         }
#     } elsif ( $check_for eq 'NOT_EXISTS' ) {
#         foreach my $file (@$files) {
#             if ( -e $file ) {
#                 warn "$file already exists\n";
#                 $not_exists = 0;
#             }
#         }
#     } else {
#         die "Subroutine 'check_files_exist' requires argument 'EXISTS' or 'NOT_EXISTS'.\n";
#     }

#     if ( $check_for eq 'EXISTS' ) {
#         return $all_exist;
#     } else {
#         return $not_exists;
#     }
# }


