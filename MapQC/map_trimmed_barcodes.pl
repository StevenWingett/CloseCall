#!/usr/bin/perl

##############################################################
#Takes a FASTQ file in which the first 26bps constitute
#a barcode. These are converted into a pseudo reference
#genome and a FASTA file
#
#The script then builds the indices and maps the barcode against these
#However, the barcode are trimmed by 3bp at either end before mapping
##############################################################

use strict;
use warnings;

use Data::Dumper;

unless(@ARGV){
	die "Please specify a file to process.\n";
}

foreach my $file (@ARGV){

	print "Processing $file\n";
	 
	if ( $file =~ /\.gz$/ ) {
		open( IN, "zcat $file |" ) or die "Couldn't read file '$file' : $!";
	} else {
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}
	
	my $genome_outfile = $file;
	$genome_outfile =~ s/\.barcode_pass\.fastq\.gz$//;
	$genome_outfile .= '.barcode_genome.fa';
	
	my $fasta_outfile = $file;
	$fasta_outfile =~ s/\.barcode_pass\.fastq\.gz$//;
	$fasta_outfile .= '.barcode_list.fasta';
	
	open(GENOME, '>', $genome_outfile) or die "Could write to '$genome_outfile' : $!";
	open(FASTA, '>', $fasta_outfile) or die "Could write to '$fasta_outfile' : $!";	

	print GENOME '>' . "$file\n";    #FASTA header
	
	my $n_string = 'N' x 24;
	
	my $i = 0;
	my %unique_barcodes;
	while(<IN>){
		my $barcode = substr(scalar <IN>, 0, 26);
		(undef) = scalar <IN>;
		(undef) = scalar <IN>;
		$i++;
		
		next if($barcode =~ /N/);  #Ignore reads containing NS

		
		
		next if(exists $unique_barcodes{$barcode});
		
		$unique_barcodes{$barcode} = '';
				
		print GENOME $barcode . $n_string . "\n";
		
		my $barcode_trimmed = substr($barcode, 3, 20);    #Trim the barcode prior to mapping
		print FASTA '>' . "$i.$barcode\n" . "$barcode_trimmed\n";

	}
	
	close IN or die "Could not close '$file'\n";
	close GENOME or die "Could not close '$genome_outfile' : $!";
	close FASTA or die "Could not close '$fasta_outfile' : $!";
	
	#Build Bowtie2 indices
	#my $build_command_cluster_job_id = generateRandomString(6);
	#$build_command_cluster_job_id = "$file.build.$build_command_cluster_job_id";
	
	my $build_command = "bowtie2-build $genome_outfile $genome_outfile.bowtie2_indices >/dev/null";
	#my $build_command_cluster = "qsub -l h_vmem=10G -pe cores 4 -o $build_command_cluster_job_id.out -N $build_command_cluster_job_id $build_command";
	
	#print "Cluster Job Build Index: '$build_command_cluster'\n";
	print "Build Index: '$build_command'\n";
	!system($build_command) or die "Could not run '$build_command'\n";
	
	
	#Map using Bowtie2
	#my $bowtie2_command_cluster_job_id = generateRandomString(6);
	#$bowtie2_command_cluster_job_id = "$file.bowtie2.$bowtie2_command_cluster_job_id";
	
	my $bowtie2_command = "bowtie2 -f -a -x $genome_outfile.bowtie2_indices $fasta_outfile -S $fasta_outfile.sam";
	#my $bowtie2_command_cluster = "qsub -l h_vmem=10G -pe cores 4 -o $bowtie2_command_cluster_job_id.out -N $bowtie2_command_cluster_job_id -hold_jid $build_command_cluster_job_id $bowtie2_command";
	
	#print "Cluster Job Map Data: '$bowtie2_command_cluster'\n";
	#!system($bowtie2_command_cluster) or die "Could not run '$bowtie2_command_cluster'\n";
	
	#print "Cluster Job Map Data: '$bowtie2_command_cluster'\n";
	#!system($bowtie2_command_cluster) or die "Could not run '$bowtie2_command_cluster'\n";
	
	print "Map Data: '$bowtie2_command'\n";
	!system($bowtie2_command) or die "Could not run '$bowtie2_command'\n";
		
}


print "Processing complete.\n";


exit (0);


###################################################################
#Subroutines
###################################################################


#Create a random letter string
#Takes a number for the length of
#the string and returns the string
sub generateRandomString{
	my $length = $_[0];
	
	my @chars = ("A".."Z", "a".."z");
	my $string;
	$string .= $chars[rand @chars] for 1..$length;
	
	return $string;
}
