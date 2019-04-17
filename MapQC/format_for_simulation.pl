#!/usr/bin/perl

##########################################################################
#Perl script to tidy up a 'uniques' file generated by Anaconda to prepare
#it for the Monte Carlo simulation (and possibly other analysis)
##########################################################################

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/../";
use File::Basename;
use Getopt::Long;
use module_anaconda;

use Data::Dumper;


#Get user input
#Option variables
my %config = (repeats => '');
my $config_result = GetOptions(
					"repeats=s" => \$config{repeats}, 
			      );
die "Could not parse options" unless ($config_result);

die "Specify files to process.\n" unless(@ARGV);
die "Specify the --repeats file.\n" unless(hasval($config{repeats}));


#Read repeats file to identify repeats positions
print "Reading repeats file '$config{repeats}\n";
if ( $config{repeats} =~ /.*\.gz$/ ) {
	open( REPEATS, "zcat $config{repeats} |" ) or die "Cannot open '$config{repeats}' : $!";
} else {
	open( REPEATS, $config{repeats} ) or die "Cannot open '$config{repeats}' : $!";
}

my %summary_counter = (Reads_Processed => 0, Reads_Retained => 0, Percent_Reads_Retained => 0);

scalar <REPEATS>;    #Ignore header 

my %repeat_list;  # %{"Csome\tStart\End"} = ''  
my %first_repeat_sequence;    # %{id} = Csome\tStart\tEnd of 1 representative repeat, closest to the "start of the genomic sequence"
while (<REPEATS>) {
	my $line = $_;
	chomp $line;
	next if ($line =~ /^\s*$/);

	my ($csome, $start, $end, $strand, $repeat_name) = split(/\t/, $line);		
	my $cse = "$csome\t$start\t$end";

	$repeat_list{$cse} = '';

	if( exists $first_repeat_sequence{$repeat_name} ) {   #Create %first_repeat_sequence
		if( ($first_repeat_sequence{$repeat_name} cmp $cse) == 1 ){
			$first_repeat_sequence{$repeat_name} = $cse;
		}
	}else{
		$first_repeat_sequence{$repeat_name} = $cse;
	}
}

close REPEATS or die "Could not close '$config{repeats}' : $!";

#print Dumper \%first_repeat_sequence;
#print Dumper \%repeat_list;


#Process file(s)
print "Processing files:\n";
my @files = deduplicate_array(@ARGV);
foreach my $file (@files){
	print "\t$file\n";

	if( $file =~ /\.gz$/ ) {
    	open( INTERACTIONS, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
    	open( INTERACTIONS, '<', $file ) or die "Could not read '$file' : $!";
	}

	my %deduplicated_results;
	my $header = scalar <INTERACTIONS>;


	#Convert repeats to the position the repeat first occurs in the genome
	while(<INTERACTIONS>){
		my $line = $_;
		chomp $line;
		next if ($line =~ /^\s*$/);
		$summary_counter{Reads_Processed}++;

		my @line_elements = split(/\t/, $line);
		my $barcode_id = $line_elements[1];
		my $csome = $line_elements[2];
		my $start = $line_elements[3];
		my $end = $line_elements[4];
		my $name = $line_elements[7];
		my @name_elements = split(/\./, $name);    #Remove numerical indices
		my $first_index_id = shift @name_elements;   #Use later
		shift @name_elements;
		$name = join('.', @name_elements);

		my $lookup_id = "$csome\t$start\t$end";

		my $is_repeat_flag = 0;
		if( exists $repeat_list{$lookup_id} ){    #This is a repeat
			if( exists $first_repeat_sequence{$name} ){
				
				#Don't change the position any more			
				#my $new_position = $first_repeat_sequence{$name};
				#my ($new_csome, $new_start, $new_end) = split(/\t/, $new_position);
				#$line_elements[2] = $new_csome;
				#$line_elements[3] = $new_start;
				#$line_elements[4] = $new_end;
				
				#Remove numerical index for repeats
				$line_elements[7] = join('.', @name_elements);
				$is_repeat_flag = 1;
			}else{
				die "Could not find '$name' in first_repeat_sequence hash.\n";    #This should not happen!
			}
		}else{
		
			#Format the id so the different exons are no longer recorded 
			#e.g. 1.2.GeneA -> 1.GeneA
			$line_elements[7] = join('.', @name_elements);
			$line_elements[7] = "$first_index_id.$line_elements[7]";  
		}

 	
		#my $results_string = join("\t", @line_elements);


		#De-duplicate the dataset (thereby preventing inter-gene interactions)
		#The same feature within the same barcode group is considered a duplicate 
		if( (exists $deduplicated_results{"$barcode_id\t$name"}) and ($is_repeat_flag == 0) ){    #Duplicate but not repeat
			my @result_in_hash = split(/\t/, $deduplicated_results{"$barcode_id\t$name"});

			#Edit the start and end position so they span the maximum distance
			#covered by the duplicates
			if($result_in_hash[3] < $line_elements[3]){   #Check start position
				$line_elements[3] = $result_in_hash[3]
			} 

			if($result_in_hash[4] > $line_elements[4]){   #Check end position
				$line_elements[4] = $result_in_hash[4]
			} 
			$deduplicated_results{"$barcode_id\t$name"} = join("\t", @line_elements);

		}else{    #Fist time observed in dataset, or repeat
			$deduplicated_results{"$barcode_id\t$name"} = join("\t", @line_elements);
		}

	}
	close INTERACTIONS or die "Could not close filehandle on '$file' : $!";

	#print Dumper \%deduplicated_results;

	#Print out the results
	my %sorted_output_results;    #Sort the output by read ID %{read_id} = result for printing
	foreach my $result (sort keys %deduplicated_results){
		my $result = $deduplicated_results{$result};
		my $read_id = (split(/\t/, $result))[0];
		$sorted_output_results{$read_id} = $result;
	}


	my $outfile = $file;
	$outfile =~ s/\.uniques\.txt\.gz//;
	$outfile .= ".simulation_formatted.txt.gz";

	open(OUT, "| gzip -c - >  $outfile") or die "Could not write to '$outfile' : $!";
	print OUT $header;
	foreach my $read_id (sort {$a <=> $b} keys %sorted_output_results){
		print OUT "$sorted_output_results{$read_id}\n";
		$summary_counter{Reads_Retained}++;
	}
	close OUT or die "Could not close filehandle on '$outfile' : $!";

	#Print out the summary results file
	$summary_counter{Percent_Reads_Retained} = calc_perc($summary_counter{Reads_Retained}, $summary_counter{Reads_Processed}, 2);
	my $summary_file = "$outfile.summary.txt";
	open (SUMMARY, '>', $summary_file) or die "Could not write to summary file '$summary_file' : $!";
	print SUMMARY "File\tProcessed\tRetained\tPercent_Retained\n";
	print SUMMARY "$file\t$summary_counter{Reads_Processed}\t$summary_counter{Reads_Retained}\t$summary_counter{Percent_Reads_Retained}\n";
	close SUMMARY or die "Could not close filehandle on '$summary_file' : $!";

}

print "Edited files for simulation.\n";

exit (0);
