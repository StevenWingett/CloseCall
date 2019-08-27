#!/usr/bin/perl

################################################################
#A Perl script to collate the results of multiple simulation files
################################################################

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../";
use Data::Dumper;


#Option variables
my %config = (
	observed => '',
	prefix => '',
	simulations => '',
	sleep => 0,
	help => ''
);
	
my $config_result = GetOptions(
	"prefix=s" => \$config{prefix},
	"simulations=i"	=> \$config{simulations},
	"sleep=i" => \$config{sleep},
    "help"     	=> \$config{help},
);
die "Could not parse options\n" unless ($config_result);

if($config{sleep}){    #Ensure input files all generated before processing
	sleep ($config{sleep});
}

unless(@ARGV){
	die "Please specify a file to process.\n";
}

my @files = dedup_array_keep_order(@ARGV);


#Read in data
my $total_number_simulations = 0;
my %collated_results;    #%{interaction} = @(obs_freq, sim_Avg_Freq, sim_score);
my %previous_collated_results;  #As collated results but 1 file 'behind'

print "Collating results:\n";

my $obs_sim_changes_outfile = "obs_sim_changes.txt";
my $p_val_changes_outfile = "p_value_changes.txt";
if( $config{prefix} ne ''){
	$obs_sim_changes_outfile = "$config{prefix}.$obs_sim_changes_outfile";
	$p_val_changes_outfile = "$config{prefix}.$p_val_changes_outfile";
}


open(OBS_SIM, '>', $obs_sim_changes_outfile) or die "Could not write to '$obs_sim_changes_outfile' : $!";
open(PVAL, '>', $p_val_changes_outfile) or die "Could not write to '$p_val_changes_outfile' : $!";

my $files_processed = 1;
foreach my $file (@files){
	print "\treading in '$file'\n";

	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	}else{
		open(IN, '<', $file) or die "Could not open '$file' : $!";
	}
	
	my $number_simulations;
	scalar <IN>;    #Ignore header
	while(<IN>){
		my $line = $_;
		chomp $line;
		my ($feature1, $feature2, $obs_Freq, $sim_Avg_Freq, $obs_Sim, $sim_score, $p_value) = split(/\t/, $line);
		my $interaction = "$feature1\t$feature2";


		#Determine number of simulations
		if( (!defined $number_simulations) and ($sim_score != 0) ){
			$number_simulations = sprintf("%.0f", $sim_score / $p_value);     #Convert to integer to prevent later precision problems
		}

		if(exists $collated_results{$interaction}){
			$ { $collated_results{$interaction} }[1] = ${ $collated_results{$interaction} }[1] + $sim_Avg_Freq;
			$ { $collated_results{$interaction} }[2] = ${ $collated_results{$interaction} }[2] + $sim_score;
		}else{
			$ { $collated_results{$interaction} }[0] = $obs_Freq;
			$ { $collated_results{$interaction} }[1] = $sim_Avg_Freq;
			$ { $collated_results{$interaction} }[2] = $sim_score;
		} 
	}
	close IN or die "Could not close $file : $!";

	if( !defined($number_simulations)){
		if( $config{simulations} ne '' ){
			die "Could not determine the number of simulations performed.\n";
		}
		print "Setting the number of simulations PER FILE to $config{simulations}\n";
		$number_simulations = $config{simulations};
	}
	
	if(defined $number_simulations and ( $config{simulations} ne '') ){
			
		unless($number_simulations == $config{simulations}){    #DOES == always work here??!
			warn "When reading '$file':\n";		
			warn "specified --simulations per input file ('$config{simulations}') does not correspond to that calculated ('$number_simulations')!\n";			
		}
	}	
	$total_number_simulations = $total_number_simulations + $number_simulations; 
	
	
	#Record how Observed/Simulated and P-values have changed with the addition of this data
	#print Dumper \%collated_results;
	#print Dumper \%previous_collated_results;
	
	
	chart_progress($number_simulations) if($files_processed > 1);    #Start when 2 files have been processed	
	%previous_collated_results = ();
	foreach my $interaction (keys %collated_results){
		my @array = @{ $collated_results{$interaction} };
		$previous_collated_results{$interaction} = [@array];
	}
	#$previous_collated_results{A}++;
	$files_processed++;
}

close OBS_SIM or die "Could not close filehandle on '$obs_sim_changes_outfile' : $!";
close PVAL or die "Could not close filehandle on '$p_val_changes_outfile' : $!";


#Write collated data to an output file
my $outfile = "simulation_collated_data.txt.gz";
if( $config{prefix} ne ''){
	$outfile = "$config{prefix}.$outfile";
}

open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";
print OUT "Name_Feature\tName_Feature2\tObserved_Frequency\tSimulation_Average_Frequency\tObserved/Simulation\tSimulation_Score\tP_Value\n";


print "Writing collated data to '$outfile'\n";
foreach my $interaction (sort { $collated_results{$a} <=> $collated_results{$b} } keys %collated_results) {
	my $obs_Freq = ${ $collated_results{$interaction} }[0];
	my $sim_Avg_Freq = ${ $collated_results{$interaction} }[1] / scalar(@files);
	my $sim_score  = ${ $collated_results{$interaction} }[2];

	my $obs_sim;
	if($sim_Avg_Freq == 0){
		$obs_sim = "NA";
	}else{
		$obs_sim = $obs_Freq / $sim_Avg_Freq;
	}
	my $p_value = $sim_score / $total_number_simulations;
	print OUT "$interaction\t$obs_Freq\t$sim_Avg_Freq\t$obs_sim\t$sim_score\t$p_value\n";
}

print "Creating collation summary graphs\n";

my $command = "Rscript $Bin/RScripts/monteCarloQC.r >/dev/null";
!system($command) or die "Could not create simulation QC graphs with command: '$command'\n";

#Zip the simulation progress files now (could not be done earler since R used 'fread' to read in the data
$command = "gzip $obs_sim_changes_outfile";
!system($command) or die "Could not perform command: '$command'\n";
$command = "gzip $p_val_changes_outfile";
!system($command) or die "Could not perform command: '$command'\n";

print "Processing complete\n";

exit (0);








############################################################################
#Subroutines
############################################################################

#Subroutine chart_progress
#Uses the %previous_collated_results and %collated_results to record how
sub chart_progress{
	
	
	#print "###############################\n\n\n";
	#print Dumper \%collated_results;
	#print Dumper \%previous_collated_results;

	my $number_simulations = $_[0];
	my $obs_sim_line_to_print = '';
	my $p_value_line_to_print = '';
	my $default_min = 1 / 1_000_000;    #Prevent division by zero errors
	
	foreach my $interaction (keys %previous_collated_results){
				
		#Compare obs/sim_Avg_Freq
		my $previous_obs = ${ $previous_collated_results{$interaction} }[0];
		my $previous_sim = ${ $previous_collated_results{$interaction} }[1] / ($files_processed - 1 );
		$previous_sim = $default_min if($previous_sim == 0);	
		my $previous_obs_sim = $previous_obs / $previous_sim ;
		
		
		my $current_obs = ${ $collated_results{$interaction} }[0];
		my $current_sim = ${ $collated_results{$interaction} }[1] / $files_processed;
		$current_sim = $default_min if($current_sim == 0);
		my $current_obs_sim = $current_obs / $current_sim ;
		
		my $obs_sim_change = $current_obs_sim / $previous_obs_sim;
		$obs_sim_line_to_print .= "$obs_sim_change\t";
		
		
		#Compare P values
		my $previous_sim_score = ${ $previous_collated_results{$interaction} }[2];
		my $current_sim_score = ${ $collated_results{$interaction} }[2];
		
		
		my $previous_p_value =  $previous_sim_score / ( ($files_processed - 1) * $number_simulations);
		$previous_p_value = $default_min if($previous_p_value == 0);    #Set minimum value at limit of detection to get division by zero errors
				
		my $current_p_value = $current_sim_score / ($files_processed * $number_simulations);			
		my $p_value_change = $current_p_value / $previous_p_value;
		$p_value_line_to_print .= "$p_value_change\t";

	}
		$obs_sim_line_to_print =~ s/\t$/\n/;
		$p_value_line_to_print =~ s/\t$/\n/;
		
		print OBS_SIM $obs_sim_line_to_print;
		print PVAL $p_value_line_to_print;
}



############################
#Subroutine takes an array, returns the
#de-duplicated array but maintains the original order
#(if a duplicate is found in multiple position in an array,
#its first occurrence is used) 
sub dedup_array_keep_order {
	my %already_seen;
	my @dedup_aray;
	
	foreach my $element (@_){
		unless(exists $already_seen{$element}){
			push (@dedup_aray, $element);
			$already_seen{$element} = '';
		}
	}
	return @dedup_aray;
}








