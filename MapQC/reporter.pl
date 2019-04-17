#!/usr/bin/perl

##############################################################################
#Perl script takes a list of the Anaconda input files and processes all the
#summary files which should be in the same folder
##############################################################################

use strict;
use warnings;

use Data::Dumper;

#Check input ok
die "Please specify name of input file(s) passed to Anaconda to process summary files.\n" unless (@ARGV);

print "Generating summary files for:\n";

my @summary_extensions = qw( 	.summary_sequence_present.txt
								.present.fastq.gz.problem_barcodes_summary_file.txt
								.barcode_list.fasta.sam.group_barcodes_summary.txt
								.barcode_pass.fastq.gz.assigned_barcodes_summary.txt
								.trimmed.fastq.anaconda_mapper_summary_file.txt
								.trimmed.sam.anaconda_mapper_summary_file.txt
								.filtered_reads.bam.summary_condensed_data.txt
								.condensed_format.txt.on_target_summary.txt
								.on_target.txt.gz.uniques.summary.txt
								.simulation_formatted.txt.gz.summary.txt
						);    # i.e. Input_Filename.Summary_Extension



my %results;    #%{Input_File}->{Pipeline_Step}->{Header} = result 
#Open each summary file in turn and write the results into a hash


foreach my $file(@ARGV){
	print "\t$file\n";
	foreach my $extension(@summary_extensions){
		my $summary_file = $file . $extension;
			open(IN, '<', $summary_file) or die "Could not open '$summary_file' : $!";
			

			my $header_line = scalar <IN>;
			chomp $header_line;
			my @headers = split(/\t/, $header_line);
			my $value_line = scalar <IN>;
			chomp $value_line;
			my @values = split(/\t/, $value_line);

			for my $i (0 .. $#headers) {
				my $header = $headers[$i];
				my $value = $values[$i];
				$results{$file}->{$extension}->{$header} = $value;
			}	
			close IN or die "Could not close filehandle on '$summary_file' : $!";
	}
}

#print Dumper \%results;

#Create Anaconda summary output files
my @summary_categories = qw (    
								Sample:Sample_Name
								Fixed_Seq_Check:Processed
								Fixed_Seq_Check:Total_Matching
								Fixed_Seq_Check:Percent_Matching
								Barcodes:LC_Adapter_Pass
								Barcodes:%LC_Adapter_Pass
								Barcodes:Unique_Barcodes
								Barcodes:Barcode_Groups
								Barcodes:Reads_(C)_allocated_to_Barcode_Groups_(H)
								Mapping:Unique_Map
								Mapping:Multi-map_repeat
								Mapping:Mult-map_non_repeat
								Mapping:Unmapped
								Mapping:%Mapped
								Parsing:%PCR_Uniques
								Parsing:%Genes_Repeats
								Parsing:%Anaconda_Uniques
								Parsing:Anaconda_Reads
							);    #Summary outfile has 2 headers. Header1/Header2 separated by a colon above



#%{Input_File}->{Pipeline_Step}->{Header} = result 
foreach my $file (@ARGV){

	my %summary_results_for_printing;    # %{Header1:Header2} = result
	$summary_results_for_printing{'Sample:Sample_Name'} = $file;
	$summary_results_for_printing{'Fixed_Seq_Check:Processed'} = $results{$file}->{'.summary_sequence_present.txt'}->{'Processed'};
	$summary_results_for_printing{'Fixed_Seq_Check:Total_Matching'} = $results{$file}->{'.summary_sequence_present.txt'}->{'Total_Matching'};
	$summary_results_for_printing{'Fixed_Seq_Check:Percent_Matching'} = $results{$file}->{'.summary_sequence_present.txt'}->{'Percent_Matching'};
	$summary_results_for_printing{'Barcodes:LC_Adapter_Pass'} = $results{$file}->{'.present.fastq.gz.problem_barcodes_summary_file.txt'}->{'Pass'};
	$summary_results_for_printing{'Barcodes:%LC_Adapter_Pass'} = $results{$file}->{'.present.fastq.gz.problem_barcodes_summary_file.txt'}->{'Percent_Pass'};
	$summary_results_for_printing{'Barcodes:Unique_Barcodes'} = $results{$file}->{'.barcode_list.fasta.sam.group_barcodes_summary.txt'}->{'Unique_Barcodes'};
	$summary_results_for_printing{'Barcodes:Barcode_Groups'} = $results{$file}->{'.barcode_list.fasta.sam.group_barcodes_summary.txt'}->{'Unique_Barcode_Groups'};
	$summary_results_for_printing{'Barcodes:Reads_(C)_allocated_to_Barcode_Groups_(H)'} = $results{$file}->{'.barcode_pass.fastq.gz.assigned_barcodes_summary.txt'}->{'Reads_assigned_barcode'};					
	$summary_results_for_printing{'Mapping:Unique_Map'} = $results{$file}->{'.trimmed.sam.anaconda_mapper_summary_file.txt'}->{'Unique_map'};
	$summary_results_for_printing{'Mapping:Multi-map_repeat'} = $results{$file}->{'.trimmed.sam.anaconda_mapper_summary_file.txt'}->{'Multi_map_allowed'};
	$summary_results_for_printing{'Mapping:Mult-map_non_repeat'} = $results{$file}->{'.trimmed.sam.anaconda_mapper_summary_file.txt'}->{'Multi_map_rejected'};						
	$summary_results_for_printing{'Mapping:Unmapped'} = $results{$file}->{'.trimmed.sam.anaconda_mapper_summary_file.txt'}->{'Unmapped'};						
	$summary_results_for_printing{'Mapping:%Mapped'} = $results{$file}->{'.trimmed.sam.anaconda_mapper_summary_file.txt'}->{'PC_Total_Allowed'};						
	$summary_results_for_printing{'Parsing:%PCR_Uniques'} = $results{$file}->{'.filtered_reads.bam.summary_condensed_data.txt'}->{'PC_Unique_Reads'};																				
	$summary_results_for_printing{'Parsing:%Genes_Repeats'} = $results{$file}->{'.condensed_format.txt.on_target_summary.txt'}->{'Percent_On_Target'};																		
	$summary_results_for_printing{'Parsing:%Anaconda_Uniques'} = $results{$file}->{'.on_target.txt.gz.uniques.summary.txt'}->{'Percent_Retained'};																		
	$summary_results_for_printing{'Parsing:Anaconda_Reads'} = $results{$file}->{'.simulation_formatted.txt.gz.summary.txt'}->{'Retained'};								# 
	
	#print Dumper \%summary_results_for_printing;

	my $summary_outfile = "$file.anaconda_summary_results.txt";
	open(OUT, '>', $summary_outfile);

	my $header1 = '';
	my $header2 = '';
	foreach my $category (@summary_categories){    #Print Headers
		$header1 = $header1 . (split(/:/, $category))[0] . "\t";
		$header2 = $header2 . (split(/:/, $category))[1] . "\t";
	}
	$header1 =~ s/\t$/\n/;
	$header2 =~ s/\t$/\n/;
	print OUT $header1 . $header2;

	my $results_line = '';
	foreach my $category (@summary_categories){    #Print results
		$results_line = $results_line . $summary_results_for_printing{$category} . "\t";
	}
	$results_line =~ s/\t$/\n/;
	print OUT $results_line;

	close OUT or die "Could not close filehandle on '$summary_outfile' : $!";
}

print "Generated Anaconda summary files.\n";

exit (0);
