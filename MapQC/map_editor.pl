#!/usr/bin/perl

###############################################################
#Opens a SAM file and writes unique mapping
#reads to the output file.
###############################################################

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Math::Round;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;

############################
#Get and validate user input
#Option variables
my %config = (
    repeats => '',
);

my $config_result = GetOptions(
    "repeats=s"   => \$config{repeats},
);
die "Could not parse options.\n" unless ($config_result);

die "Please specify file(s) to process.\n" unless (@ARGV);
unless(hasval($config{repeats})){
	die "Please specify a formatted repeats file (created at the start of the Anaconda pipeline.\n";
}


######################################################
#Process the file of pre-determined repeat elements
print "Reading repeats file '$config{repeats}\n";
if ( $config{repeats} =~ /.*\.gz$/ ) {
	open( REPEATS, "zcat $config{repeats} |" ) or die "Cannot open '$config{repeats}' : $!";
} else {
	open( REPEATS, $config{repeats} ) or die "Cannot open '$config{repeats}' : $!";
}

scalar <REPEATS>;    #Ignore header 

my %repeat_sequences;    #%{csome\t10kb_region}->{start} = end
while (<REPEATS>) {

	my $csome = ( split /\t/ )[0];
    my $start = ( split /\t/ )[1];
    my $end = ( split /\t/ )[2];
	my $strand = ( split /\t/ )[3];
	add_bait($csome, $start, $end, \%repeat_sequences);
}
close REPEATS or die "Could not close '$config{repeats}' : $!";


#print Dumper \%repeat_sequences;



#######################################################
#Read in the SAM file and write to the BAM file allowed
#reads. Report the results
foreach my $file (@ARGV){

	print "Processing mapped data '$file'\n";
	
	#Prepare summary file and counter
	my $summary_filename = "$file.anaconda_mapper_summary_file.txt";
	open(SUMMARY, '>', $summary_filename) or die "Could not open '$summary_filename' : $!";
	my @summary_categories = qw( Reads_Processed Unique_map Multi_map_allowed 
		Multi_map_rejected Unmapped  PC_Unique_map PC_Multi_map_allowed 
		PC_Multi_map_rejected PC_Unmapped Total_Allowed Total_Rejected 
		PC_Total_Allowed);
	my %summary_counter;
	foreach my $category ( @summary_categories ){
		$summary_counter{$category} = 0;
	}
	
	#Open input file and write to output files
	if( $file =~ /\.bam$/ ) {
		open( IN, "samtools view -h $file |" ) or die "Couldn't read '$file' : $!";
	}elsif( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }
	
	my $outfile = $file;
	$outfile =~ s/\.trimmed\.sam$//;
	$outfile .= '.filtered_reads.bam';

	open(OUT, "| samtools view -bSh 2>/dev/null - > $outfile"  ) or die "Could not write to '$outfile' : $!";


	my %id_category;    # %{ID} = category (needed since HISAT2 may produce multiple alignment lines for a single read)
	while(<IN>){
		my $read = $_;
		chomp $read; 
		
		if( substr($read, 0, 1) eq '@'){   #SAM headers
			print OUT "$read\n";
			next;
		}
		
		my ($id) = split(/\t/, $read);
		$id = (split(/_/, $id))[1];    #Remove barcode prefix

		#If the read id is already valid, do not process again
		#i.e. this read is already represented in the ouput file
		#It should not be written there again - the same read aligned multiple times
		if( exists ($id_category{$id}) ){
			if( ($id_category{$id} eq 'Unique_map') or ($id_category{$id} eq 'Multi_map_allowed') ){
				next;
			}
		}

		
		#Did the read map?
		my $sam_flag = (split(/\t/, $read))[1];
		if($sam_flag & 0x4){
			$id_category{$id} = 'Unmapped';
			next;
		}
					
		#No secondary alignments
		if($sam_flag & 0x100){    #Secondary alignment flag
			$id_category{$id} = 'Multi_map_rejected';
			next;
		}
		
		#Rescue alignments with scores better than all others - i.e.
		#consider these unique alignments
		my $current_alignment_score = (split(/\t/, $read))[11];
		my $best_other_alignment_score = (split(/\t/, $read))[12];
		
		unless( $current_alignment_score =~ /^AS:i:-*\d+/ ){
			die "Alignment score in read did not contain proper alignment score format 'AS:i:' :\n $read\n";	
		}
		$current_alignment_score = substr($current_alignment_score, 5);
			
		unless( $best_other_alignment_score =~ /^ZS:i:-*\d+/){    #This is a unique alignment (sometimes ZS:i: is missing, sometimes there is no digit after the i)
				$id_category{$id} = 'Unique_map';
				print OUT "$read\n";
				next;
		}
		$best_other_alignment_score = substr($best_other_alignment_score, 5);
		
		if($current_alignment_score > $best_other_alignment_score){
		#	$id_category{$id} = 'Unique_map';
		#	print OUT "$read\n";
			$id_category{$id} = 'Multi_map_rejected';    #New more stringent mapping parameters
			next;
		}else{	
			$id_category{$id} = 'Multi_map_rejected';
			#No 'next' here - check below whether this maps to a pre-defined region
		}
				
		#Primary alignments that map to pre-specified regions are allowed
		my ($csome, $midpoint) = samMidPoint($read);
		my $mapped_to_repeat_flag = coord2bin($csome, $midpoint, \%repeat_sequences);
		$mapped_to_repeat_flag = 1 unless($mapped_to_repeat_flag eq '0');
		if($mapped_to_repeat_flag){   #Mapped to repeat
			print OUT "$read\n";
			$id_category{$id} = 'Multi_map_allowed';	
		}
	}
		
	close OUT or warn "Could not close filehandle on '$outfile' : $!"; 

	#Print summary results
	foreach my $id (keys %id_category){
		my $category = $id_category{$id};
		$summary_counter{$category}++;
		$summary_counter{Reads_Processed}++;
	}

	$summary_counter{PC_Unique_map} = calc_perc($summary_counter{Unique_map}, $summary_counter{Reads_Processed}, 2);
	$summary_counter{PC_Multi_map_allowed} = calc_perc($summary_counter{Multi_map_allowed}, $summary_counter{Reads_Processed}, 2);
	$summary_counter{PC_Multi_map_rejected} = calc_perc($summary_counter{Multi_map_rejected}, $summary_counter{Reads_Processed}, 2);
	$summary_counter{PC_Unmapped} = calc_perc($summary_counter{Unmapped}, $summary_counter{Reads_Processed}, 2);
	$summary_counter{Total_Allowed} = $summary_counter{Unique_map} + $summary_counter{Multi_map_allowed};
	$summary_counter{Total_Rejected} =$summary_counter{Multi_map_rejected} + $summary_counter{Unmapped};
	$summary_counter{PC_Total_Allowed} = calc_perc($summary_counter{Total_Allowed}, $summary_counter{Reads_Processed}, 2);
	
	print SUMMARY "File";
	foreach my $category (@summary_categories){
		print SUMMARY "\t$category";
	}
	print SUMMARY "\n";
	print SUMMARY "$file";
	foreach my $category (@summary_categories){
		print SUMMARY "\t$summary_counter{$category}";
	}
	print SUMMARY "\n";
	close SUMMARY or die "Could not close '$summary_filename' : $!";
	
	#print Dumper \%summary_counter;
	
}


print "Processing complete.\n";

exit (0);
