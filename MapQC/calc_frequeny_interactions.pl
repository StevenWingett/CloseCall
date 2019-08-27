#!/usr/bin/perl

################################################################
#A Perl script that takes an Anaconda file, produces summary
#statistics and reports the observed/expected value for each
#gene interaction
#
#Anaconda file input format:
#Read_ID	
#barcode_id Chromosome	
#Feature_Start	
#Feature_End	
#Feature_Strand	
#Feature_ID	
#Feature_Name
################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;


use Data::Dumper;


###########################################################################
#Verify user input

#Option variables
my %config = (
	help => ''
);
	
my $config_result = GetOptions(
    "help"     	=> \$config{help},
);
die "Could not parse options\n" unless ($config_result);

if ( $config{help} ) {    #Print help and exit
    print while (<DATA>);
    exit(0);
}

unless(@ARGV){
	warn "Please specify file(s) to process.\n";
	print while (<DATA>);
    die "\n";
}

my @files = deduplicate_array(@ARGV);
if(scalar @files > 1){
	warn "You have specified multiple files to process which may take some time (maybe submit in parrallel?)\n";
}

print "Determining frequency of interactions\n";


##########################################################################
#Read in data
my %barcodes_genes;    #%{barcode_ID} = @Genes
foreach my $file (@files){

	print "Reading in file '$file'\n";

	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	}else{
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}
	
	scalar <IN>;    #Ignore header

	while(<IN>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		shift @line_elements;
		my $barcode_id = shift @line_elements;

		if($line_elements[4] eq 'NA'){    #Repeats have value Feature_ID = NA in Anaconda file. Make this now the same as the name for repeats - this has changed
			$line_elements[4] = $line_elements[5]; 
		}

		my $remainder = join("\t", @line_elements);

		#my (undef, $barcode_id, $csome, $start, $end) = split(/\t/, $line);
		push ( @{ $barcodes_genes{$barcode_id} }, $remainder );
	}
	      
	close IN or die "Could not close $file; : $!";
	  
	#print Dumper \%barcodes_genes;
	



	##########################################################################
	#Analyse sample data

	print "Calculating interaction frequencies\n";

	#1) Determine barcode group sizes
	calc_barcode_group_size(\%barcodes_genes, $file);

	#2) Frequency distribution of interactions per gene
	calc_freq_interactions_every_gene(\%barcodes_genes, $file);

	#3) Frequency distribution of interactions between specific genes
	calc_freq_specific_interactions(\%barcodes_genes, $file);


	#4) Determine the number total number of interacting components, and
	#then number of times they are interaction (the same gene may interact)
	#multiple times. This allows a basic obsereved/expected measure described
	#below
	my %geneID_count = calc_gene_frequency(\%barcodes_genes);
	my $total_counts = 0;
	foreach my $id (keys %geneID_count){
		$total_counts += $geneID_count{$id};
	}
	#print $total_counts;
	#print Dumper \%geneID_count;


	#######################################################################
	#Write gene-gene interactions to an output file
	#Two values for Observed/Expected are returned:
	#i) Calculated using the proportion of times a gene is observed in the dataset
	#ii) Calculated using the proportion of times a gene is observed in a virtual di-tag
	my %sample_interactions_counter = calc_freq_specific_interactions(\%barcodes_genes, $file, 1);
	
	#Determine the number of times a gene is observed in a virtual di-tag
	my $virtual_ditag_count = 0;
	my %genes_appearing_in_virtual_ditags_counter;    #%{gene_id} = count
	foreach my $interaction (keys %sample_interactions_counter){
		$virtual_ditag_count += $sample_interactions_counter{$interaction};
		my $gene1_id = (split(/\t/, $interaction))[4];
		my $gene2_id = (split(/\t/, $interaction))[10];
		$genes_appearing_in_virtual_ditags_counter{$gene1_id}++;
		$genes_appearing_in_virtual_ditags_counter{$gene2_id}++;
	}

	my $outfile = $file;
	$outfile =~ s/\.simulation_formatted\.txt\.gz$//;
	$outfile .= '.partially_compressed.txt.gz';
	open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
	print OUT "Gene1_Chromosome\tGene1_Start\tGene1_End\tGene1_Strand\tGene1_ID\tGene1_Name\t";
	print OUT "Gene2_Chromosome\tGene2_Start\tGene2_End\tGene2_Strand\tGene2_ID\tGene2_Name\t";
	print OUT "Frequency\tObserved/Expected\tVirtual_Ditags_Observed/Expected\n";

	foreach my $interaction (sort { $sample_interactions_counter{$a} <=> $sample_interactions_counter{$b}  } keys %sample_interactions_counter ){
		print OUT "$interaction\t$sample_interactions_counter{$interaction}\t";
		
		my $gene1_id = (split(/\t/, $interaction))[4];
		my $gene2_id = (split(/\t/, $interaction))[10];
		
		#i) Obs/Exp calculated using the proportion of times a gene is observed in the dataset
		my $exp = ( ($geneID_count{$gene1_id} / $total_counts) * ($geneID_count{$gene2_id} / $total_counts) ) * $virtual_ditag_count;
		my $obs_exp = $sample_interactions_counter{$interaction} / $exp;
		print OUT "$obs_exp\t";		
		
		#print "$interaction\t$sample_interactions_counter{$interaction}\n";
		#print "Gene1_Count:$geneID_count{$gene1_id}\nGene2_Count:$geneID_count{$gene2_id}\nTotal_Count:$total_counts\nVirtual_Ditag_Count:$virtual_ditag_count\n";
		#print "$exp\t$obs_exp\n\n";
			

		#ii) Obs/Exp calculated using the proportion of times a gene is oberved in a virtual di-tag
		my $vd_exp = ( ($genes_appearing_in_virtual_ditags_counter{$gene1_id} / ($virtual_ditag_count * 2)) * ($genes_appearing_in_virtual_ditags_counter{$gene2_id} / ($virtual_ditag_count * 2)) ) * $virtual_ditag_count;
		my $vd_obs_exp = $sample_interactions_counter{$interaction} / $vd_exp;
		print OUT "$vd_obs_exp\n";
	} 
	close OUT or die "Could not close filehandle on '$outfile' : $!"; 
	
	#Now process the Anaconda output file to compress features.
	#The code below simply processes the orignal anaconda file. 
	#This is not an elegant solution, but it is simple and allows us to keep the stats when
	#different repeats are counted as separate repeats
	my $anaconda_file1 = $outfile;		
	if($anaconda_file1 =~ /\.gz$/){
		open (ANACONDA_1, "zcat $anaconda_file1 |") or die "Could not read file '$anaconda_file1' : $!";
	}else{
		open(ANACONDA_1, '<', $anaconda_file1) or die "Could not open '$anaconda_file1' : $!";
	}
	scalar <ANACONDA_1>;    #Ignore headers
	
	my %results;    # %{interaction} = @(Observed Expected VD_Expected Chr1 Start1 End1 Chr2 Start2 End2)	
	while(<ANACONDA_1>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		
		#Important, possibly swap the order of elements in the line so feature1, feature2 is always 
		#the same	
		if( ($line_elements[5] cmp $line_elements[11]) == 1 ){    #Order output
			my @part1 = splice(@line_elements, 0, 6);
			my @part2 = splice(@line_elements, 0, 6);
			my @part3 = splice(@line_elements, 0, 3);
			@line_elements = ();    #Should be empty anyway
			@line_elements = (@part2, @part1, @part3);
		}
		
		my $csome1 = $line_elements[0];
		my $start1 = $line_elements[1];
		my $end1 = $line_elements[2];
		my $strand1 = $line_elements[3];
				
		my $csome2 = $line_elements[6];
		my $start2 = $line_elements[7];
		my $end2 = $line_elements[8];
		my $strand2 = $line_elements[9];
			
		my $feature1_name = $line_elements[5];
		my $feature2_name = $line_elements[11];
		my $obs = $line_elements[12];
		my $obs_exp = $line_elements[13];
		my $vd_obs_exp = $line_elements[14];
		
		#If this is a repeat region, make the chromosome the repeat name and the position 1 -> 1,000,000 
		unless($feature1_name =~ /^\d+\..+$/){
			$csome1 = $feature1_name;
			$start1 = 1;
			$end1 = 1_000_000;	
			$strand1 = 'B';
		}
		
		unless($feature2_name =~ /^\d+\..+$/){
			$csome2 = $feature2_name;
			$start2 = 1;
			$end2 = 1_000_000;
			$strand2 = 'B';			
		}
			
		my $interaction = "$feature1_name\t$feature2_name";     #Ordering output should not be necessary now
	
	
		unless ( exists $results{$interaction} ){    #Initialise
			push( @{ $results{$interaction} }, ( 0, 0, 0, $csome1, $start1, $end1, $strand1, $csome2, $start2, $end2, $strand2 ) );	
		}
			
		
		#Set the feature length so it covers the maximum observed distance
		my $csome1_previous = $results{$interaction}->[3];
		my $start1_previous = $results{$interaction}->[4];
		my $end1_previous = $results{$interaction}->[5];
		my $strand1_previous = $results{$interaction}->[6];
		
		if( ($csome1_previous eq $csome1) and ($strand1_previous eq $strand1 ) ){
			if($start1 < $start1_previous){
				$results{$interaction}->[4] = $start1
			}
			if($end1 > $end1_previous){
				$results{$interaction}->[5] = $end1
			}	
		}else{
			warn "For interaction '$interaction' either '$csome1_previous' does not match '$csome1'\n";    #This should not happen	
			warn "or '$strand1_previous' does not match '$strand1'\n";
		}
		
		my $csome2_previous = $results{$interaction}->[7];
		my $start2_previous = $results{$interaction}->[8];
		my $end2_previous = $results{$interaction}->[9];
		my $strand2_previous = $results{$interaction}->[10];
		
		if( ($csome2_previous eq $csome2) and ($strand2_previous eq $strand2 ) ){
			if($start2 < $start2_previous){
				$results{$interaction}->[8] = $start2
			}
			if($end2 > $end2_previous){
				$results{$interaction}->[9] = $end2;
			}	
		}else{
			warn "For interaction '$interaction' either '$csome2_previous' does not match '$csome2'\n";    #This should not happen	
			warn "or '$strand2_previous' does not match '$strand2'\n";
		}
		
		#Sort out the observe / expected ratios
		my $expected = $obs / $obs_exp;
		$results{$interaction}->[1] += $expected;
		my $vd_expected = $obs / $vd_obs_exp;
		$results{$interaction}->[2] += $vd_expected; 
		$results{$interaction}->[0] += $obs;
				
	}	
	close ANACONDA_1 or die "Could not close filehandle on '$anaconda_file1' : $!";
	
	my $anaconda_file2 = $anaconda_file1;
	$anaconda_file2 =~ s/\.partially_compressed\.txt\.gz$/.anaconda.txt.gz/;
	open(ANACONDA_2, "| gzip -c - > $anaconda_file2") or die "Could not write to '$anaconda_file2' : $!";
	print ANACONDA_2 "Feature1_Chromsome\tFeature1_Start\tFeature1_End\tFeature1_Strand\tFeature1_Name\t";
	print ANACONDA_2 "Feature2_Chromsome\tFeature2_Start\tFeature2_End\tFeature2_Strand\tFeature2_Name\t";
    print ANACONDA_2  "Frequency\tObserved/Expected\tVirtual_Ditags_Observed/Expected\n";
	foreach my $interaction (sort keys %results){
	
		my ($feature1_name, $feature2_name) = split(/\t/, $interaction);
	
		my $csome1 = $results{$interaction}->[3];
		my $start1 = $results{$interaction}->[4];
		my $end1 = $results{$interaction}->[5];
		my $strand1 = $results{$interaction}->[6];
		
		my $csome2 = $results{$interaction}->[7];
		my $start2 = $results{$interaction}->[8];
		my $end2 = $results{$interaction}->[9];
		my $strand2 = $results{$interaction}->[10];
	
		my $observed = $results{$interaction}->[0];
		my $expected = $results{$interaction}->[1];
		my $vd_expected = $results{$interaction}->[2];
		my $obs_exp = $observed / $expected;
		my $vd_obs_exp = $observed / $vd_expected;
			
		print ANACONDA_2 "$csome1\t$start1\t$end1\t$strand1\t$feature1_name\t";
		print ANACONDA_2 "$csome2\t$start2\t$end2\t$strand2\t$feature2_name\t";
		print ANACONDA_2 "$observed\t$obs_exp\t$vd_obs_exp\n";
	}
	close ANACONDA_2 or die "Could not close filehandle on '$anaconda_file2' : $!";
}

print "Processing complete\n";

exit (0);



##########################################################################
#Subroutines
##########################################################################


##########################################
#Subroutine: calc_gene_frequency
#Takes a barcodes_genes hash reference returns
#a hash of %{gene_id} = count
sub calc_gene_frequency{
	my ($barcodes_genes_ref) = @_;
	my %geneID_count;
		
	#Calculate freqency of occurrences of all genes
	foreach my $barcode_id (keys %{ $barcodes_genes_ref } ){
		foreach my $gene_string ( @{ $$barcodes_genes_ref{$barcode_id} } ){
			my $gene_id = (split(/\t/, $gene_string))[4];
			$geneID_count{$gene_id}++;
		}
	}
	return %geneID_count;
}


##########################################
#Subroutine: calc_freq_specific_interactions
#Takes a barcodes_genes hash reference and the name
#of the file being processed and calculates 
#frequency distribution of gene-gene interactions
#e.g:
#Occurrences    Frequency
#1				3405
#2				234
#3				7
#
#Will also accept a third parameter, which if true will cause the
#subroutine to return the tally of the specific interactions only
sub calc_freq_specific_interactions{
	my ($barcodes_genes_ref, $file, $return_tally_only) = @_;
	my %gene_gene_interations_counter;
		
	#Calculate freqency of occurrences of all genes
	foreach my $barcode_id (keys %{ $barcodes_genes_ref } ){
		my @gene_gene_interations = createDitags ( @{ $$barcodes_genes_ref{$barcode_id} } );
		
		while(@gene_gene_interations){
			my $read1 = shift @gene_gene_interations;
			my $read2 = shift @gene_gene_interations;
			$gene_gene_interations_counter{"$read1\t$read2"}++;
		}
	}

	if($return_tally_only){
		return %gene_gene_interations_counter
	}

	#Calculate frequency of these frequencies
	my %occurs_freq;
	foreach my $gene_gene_interaction (keys %gene_gene_interations_counter){
		my $occurs = $gene_gene_interations_counter{$gene_gene_interaction};
		$occurs_freq{$occurs}++;
	}

	#Print out results
	my $outfile = "$file.interactions_frequencies.txt";
	open(OUT, '>', $outfile ) or die "Could not write to '$outfile ' : $!";
	print OUT "Interaction_Occurrences\tFrequency\n";
	foreach my $occurs (sort {$a <=> $b} keys %occurs_freq){
		print OUT "$occurs\t$occurs_freq{$occurs}\n";
	}
	close OUT or die "Could not close $outfile : $!";

}


##########################################
#Subroutine: calc_freq_interactions_every_gene
#Takes a barcodes_genes hash reference and the name
#of the file being processed and calculates 
#frequency distribution of interactions per gene
sub calc_freq_interactions_every_gene{
	my ($barcodes_genes_ref, $file) = @_;
	my %gene_freq;
	
	
	#Calculate freqency of occurrences of all genes
	foreach my $barcode_id (keys %{ $barcodes_genes_ref } ){
		foreach my $gene_id ( @{ $$barcodes_genes_ref{$barcode_id} } ){
			$gene_freq{$gene_id}++;	
		}
	}
	
	#Calculate frequency of these frequencies
	my %occurs_freq;
	foreach my $gene_id (keys %gene_freq){
		my $occurs = $gene_freq{$gene_id};
		$occurs_freq{$occurs}++;
	}

	#Print out results
	my $outfile = "$file.gene_frequencies.txt";
	open(OUT, '>', $outfile ) or die "Could not write to '$outfile ' : $!";
	print OUT "Occurrences\tFrequency\n";
	foreach my $occurs (sort  {$a <=> $b} keys %occurs_freq){
		print OUT "$occurs\t$occurs_freq{$occurs}\n";
	}
	close OUT or die "Could not close $outfile : $!";
	
}


##########################################
#Subroutine: calc_barcode_group_size
#Takes a barcodes_genes hash reference and the name
#of the file being processed and calculates a barcode group
#size frequency distribution
sub calc_barcode_group_size{
	my ($barcodes_genes_ref, $file) = @_;
	my %size_freq;
	
	foreach my $barcode_id (keys %{ $barcodes_genes_ref } ){
		my $barcode_group_size =  scalar @{ $$barcodes_genes_ref{$barcode_id} };
		$size_freq{$barcode_group_size}++;
	}
	
	#Print out results
	my $outfile = "$file.barcode_group_size.txt";
	open(OUT, '>', $outfile ) or die "Could not write to '$outfile ' : $!";
	print OUT "Barcode_Group_Size\tFrequency\n";
	foreach my $barcode_group_size (sort {$a <=> $b} keys %size_freq){
		print OUT "$barcode_group_size\t$size_freq{$barcode_group_size}\n";
	}
	close OUT or die "Could not close $outfile : $!";
}



######################################
#Subroutine: createDitags
#Takes an array of reads and creates an array of the read-pair combinations
#(reports read pairs once i.e. read1-read2 but not read2-read1 as well)
sub createDitags{
	my @reads = @_;
	my @ditags;
 
	while(scalar @reads > 1){
		for (my $i=1; $i < scalar @reads; $i++){
		
			if( ($reads[0] cmp $reads[$i]) == 1 ){    #Order output
				push(@ditags, $reads[$i]);
				push(@ditags, $reads[0]);
			}else{
				push(@ditags, $reads[0]);
				push(@ditags, $reads[$i]);
			}
		}	
		shift @reads;
	}
	return @ditags;
}




__DATA__

SYNOPSIS



FUNCTION


COMMAND LINE OPTIONS


--help             Print help message and exit


Steven Wingett, Babraham Institute, Cambridge, UK (steven.wingett@babraham.ac.uk)
