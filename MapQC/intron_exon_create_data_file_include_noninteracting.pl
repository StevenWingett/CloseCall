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


##############################################################################
#Perl script takes a SAM/BAM file and identifies mapped reads, removes
#duplicates and creates a data structure of
#Barcode1	Chromosome	Position	Strand
##############################################################################

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;
use Math::Round;
use PDL::LiteF;        # loads less modules
use PDL::NiceSlice;    # preprocessor for easier pdl indexing syntax 
use PDL::Stats;


use Data::Dumper;

#Check input ok
die "Please specify files to process.\n" unless (@ARGV);
die "Please adjust configuration.\n" unless( check_files_exist(\@ARGV, 'EXISTS') );



my @summaryCategories = qw(File Reads_Processed Aligned_Reads Unique_Reads Interacting_Reads Group_Too_Large PC_Unique_Reads);

foreach my $file (@ARGV){

	my $summary_filename = "$file.summary_condensed_data.txt";
	
	open (SUMMARY, '>', $summary_filename) or die "Could not write to summary file '$summary_filename'.\n";
	print SUMMARY "File\tReads_Processed\tAligned_Reads\tUnique_Reads\tBarcode_Groups\tGroup_Too_Large_Removed_Reads\tPC_Unique_Reads\n";
	print "Processing $file\n";

	if ($file =~ /\.bam$/){
		open (IN,"samtools view -h $file |") or die "Could not read '$file' : $!";
	}else{
		open (IN, '<', $file) or die "Could not read '$file' : $!";
	}
	
	my %summaryCounter;
	foreach my $category(@summaryCategories){
		$summaryCounter{$category} = 0;
	}
	
	my %uniqueReads;  #%{barcode csome position strand} = count
	my %uniqueReads_NoBarcode;  #%{csome position strand} = count  
	while(<IN>){
		my $read = $_;
		chomp $read;

		next if substr($read, 0, 1) eq '@';    #Ignore header   
		
		$summaryCounter{Reads_Processed}++;

		my ($barcode, $samFlag) = split(/\t/, $read);
		next if($samFlag & 0x4);   #Check read aligned
		$summaryCounter{Aligned_Reads}++;

		($barcode) = split (/_/, $barcode);
		my ($csome, $pos, $strand) = sam2midcoords($read);

		$uniqueReads{"$barcode\t$csome\t$pos\t$strand"}++;
		$uniqueReads_NoBarcode{"$csome\t$pos\t$strand"}++;
	}
	
	close IN or warn "Could not close '$file' : $!";

	#print Dumper \%uniqueReads;
	
	
	#Convert the %unique reads into a new data structure:
	#%{barcode} = @(shorter read description)
	my %groupedReads;
	foreach my $identifier (keys %uniqueReads){
		my ($barcode, $csome, $pos, $strand) = split(/\t/, $identifier);
		my $shortReadIdentifier = join("\t", $barcode, $csome, $pos, $strand);    #MODIFIED THIS FOR INTRON/EXON SPLIT
		push (@{ $groupedReads{$barcode} }, $shortReadIdentifier);
	}
	
	#Determine the threshold for barcode groups considered too large
	my $max_allowed_barcode_group_size = determine_max_allowed_barcode_group_size(\%groupedReads);
	#my $max_allowed_barcode_group_size = 10000;    #Turn off threshold

	#Print the results to a data output file
	my $outfile_results = $file;
	$outfile_results =~ s/\.filtered_reads\.bam//;
	my $too_large_outfile = "$outfile_results.includes_too_large.txt";
	$outfile_results .= '.condensed_format.txt';

	open(OUT, '>', $outfile_results ) or die "Could not write to '$outfile_results' : $!";
	open(INC_TOO_LARGE, '>', $too_large_outfile ) or die "Could not write to '$too_large_outfile' : $!";
	#For the barcode group size graph (to be plotted by R script barcode_group_size.r), 
	#determine the number of bars plotted to be placed to the left of the threshold line.
	my %allowed_graph_bars;    #%{size_allowed} 
	foreach my $barcode (keys %groupedReads){
		my $barcode_group_size = scalar @{ $groupedReads{$barcode} };
		if( ($barcode_group_size > 0) and ($barcode_group_size <= $max_allowed_barcode_group_size) ){
			$allowed_graph_bars{$barcode_group_size} = '';
		}
	}
	my $number_allowed_graph_bars = scalar(keys %allowed_graph_bars);
	print INC_TOO_LARGE "Allowed_Graph_Bars\t$number_allowed_graph_bars\n";

	my $i = 1;
	foreach my $barcode (keys %groupedReads){
		my $barcode_group_size = scalar @{ $groupedReads{$barcode} };
		if($barcode_group_size > 0){     #Include non-interacting
			foreach my $shortReadIdentifier ( @{ $groupedReads{$barcode} }){    #For the R script that plots barcode group sizes
				print INC_TOO_LARGE "$shortReadIdentifier\n";    #Print an index number instead of the actual barcode
			}
			if( ($barcode_group_size <= $max_allowed_barcode_group_size) ){   #Only allowed barcode groups
				$summaryCounter{Interacting_Reads}++;
				foreach my $shortReadIdentifier ( @{ $groupedReads{$barcode} }){
					print OUT "$shortReadIdentifier\n";    #Print an index number instead of the actual barcode
				}
			}else{
				$summaryCounter{Group_Too_Large}++;
			}
		}
		$i++;
	}

	close OUT or die "Could not close filehandle on '$outfile_results' : $!";
	close INC_TOO_LARGE or die "Could not close filehandle on '$too_large_outfile' : $!";

	$summaryCounter{Unique_Reads} = keys (%uniqueReads);
	$summaryCounter{PC_Unique_Reads} = calc_perc ($summaryCounter{Unique_Reads}, $summaryCounter{Reads_Processed}, 2 );
	print SUMMARY "$file\t$summaryCounter{Reads_Processed}\t$summaryCounter{Aligned_Reads}\t$summaryCounter{Unique_Reads}\t$summaryCounter{Interacting_Reads}\t$summaryCounter{Group_Too_Large}\t$summaryCounter{PC_Unique_Reads}\n";			
	close SUMMARY or die "Could not close '$summary_filename' : $!";

	#Create data tables for a plot of read frequency including barcodes
	my %freq_tally;    #%{frequency} = tally
	foreach my $barcode_read (keys %uniqueReads){
		my $count = $uniqueReads{$barcode_read};
		$freq_tally{$count}++;
	}

	my $outfile = "$file.read_freq.txt";
	open(TALLY, '>', $outfile) or die "Could not write to '$outfile' : $!";
	print TALLY "Count_Read_Observed\tTally\n";
	foreach my $count ( sort {$a <=> $b} keys %freq_tally ){
		print TALLY "$count\t$freq_tally{$count}\n";
	}
	close TALLY or die "Could not close '$outfile' : $!";

	my $plot_command = "Rscript $Bin/R_Scripts/plot_frequency.r $outfile 'Reads Tally', 'plot.pdf' > /dev/null";
	!system($plot_command) or warn "Could not create frequency plot with command: '$plot_command'\n";

	#And again, but not including barcodes
	%freq_tally = ();    #%{frequency} = tally
	foreach my $barcode_read (keys %uniqueReads_NoBarcode){
		my $count = $uniqueReads_NoBarcode{$barcode_read};
		$freq_tally{$count}++;
	}

	$outfile = "$file.read_freq_no_barcodes.txt";
	open(TALLY, '>', $outfile) or die "Could not write to '$outfile' : $!";
	print TALLY "Count_Read_Observed\tTally\n";
	foreach my $count ( sort {$a <=> $b} keys %freq_tally ){
		print TALLY "$count\t$freq_tally{$count}\n";
	}
	close TALLY or die "Could not close '$outfile' : $!";

	$plot_command = "Rscript $Bin/R_Scripts/plot_frequency.r $outfile 'Reads Tally (no barcodes)', 'plot.pdf' > /dev/null";
	!system($plot_command) or warn "Could not create frequency plot with command: '$plot_command'\n";

}


print "Processing complete.\n";

exit (0);



#######################################################################################
#Subroutines                                                                          #
#######################################################################################


######################################################
#Subroutine: determine_max_allowed_barcode_group_size
#Takes the hash %{barcode} = @(shorter read description) and returns
#the maximum allowed barcode group size
sub determine_max_allowed_barcode_group_size{
	my $groupedReads_ref = $_[0];
	my @input;    #Array of barcode group sizes

	#print Dumper \%{ $groupedReads_ref };
	
	foreach my $barcode (  keys %{ $groupedReads_ref } ) {
		my $barcode_group_size = scalar @{ ${ $groupedReads_ref }{$barcode} };
		push (@input, $barcode_group_size);
	}

	#print Dumper \@input;
	
	my $data = pdl @input;
	my $lambda = $data->mle_poisson();
	#print "lambda: $lambda\n";

	#Determine threshold
	my $barcode_group_size = 2;
	#my $lambda = 5.22;

	while(1){    #Infinite loop may happen here?!!!!!!
		my $i = 1;
		my $group_size_factorial = 1;
		$group_size_factorial *= ++$i while $i < $barcode_group_size;    #Determine lambda factorial
	#	print "Factorial : $group_size_factorial\n";
		my $prob =  ($lambda ** $barcode_group_size) * ( (exp (-$lambda)) / $group_size_factorial );
		#print "Group size: $barcode_group_size    Prob: $prob\n";
		last if ($prob < 0.001);    #1% cutoff
		$barcode_group_size++;
	}

	my $max_allowed_barcode_group_size = $barcode_group_size - 1;
	#print "Max allowed: $max_allowed";
	
	return $max_allowed_barcode_group_size;
}





########################################
#Subroutine: sam2midcoords
#Receives a SAM format line and returns
#the middle of the read
sub sam2midcoords {
	my $read = $_[0];

	my $csome = ( split( /\t/, $read ) )[2];
	my $pos   = ( split( /\t/, $read ) )[3];
	my $seq   = ( split( /\t/, $read ) )[9];
	my $strand;

	$pos = $pos + int((length($seq) - 1) / 2);

	if ( ( ( split( /\t/, $read ) )[1] ) & 0x10 ) {    #Analyse the SAM bitwise flag to determine strand
		$strand = '-';                                 #Negative strand
	} else {
		$strand = '+';                                 #Positive strand
	}
	return ( $csome, $pos, $strand );
}

