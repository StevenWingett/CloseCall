#!/usr/bin/perl

################################################################
#A Perl script takes an Anaconda format file and assigns P-Values
#to the interactions
################################################################

use strict;
use warnings;
use Getopt::Long;
#use FindBin '$Bin';
#use lib "$Bin/../";
#use module_anaconda;
use Data::Dumper;

unless(@ARGV){
	warn "Please specify file(s) to process.\n";
    die "\n";
}

#my @files = deduplicate_array(@ARGV);
my @files = @ARGV;

my $outfile = "combined_anaconda_results.txt.gz";
open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";

print "Combining Anaconda files:\n";

my %interaction_counts;    # %{10-column_interaction_description} = @(total_observed, total_expected, total_vd_expected) 
my $header;
foreach my $file (@files){
	print "\t$file\n";

	#Determine the number of times a feature is present in the dataset
	if($file =~ /\.gz$/){
		open (ANACONDA, "zcat $file |") or die "Could not read file '$file' : $!";
	}else{
		open(ANACONDA, '<', $file) or die "Could not open '$file' : $!";
	}
	$header = scalar <ANACONDA>;    #Print header to outfile
	
	while(<ANACONDA>){
		my $line = $_;
		chomp $line;
		my @line_elements = split(/\t/, $line);
		my $vd_expected = pop @line_elements;
		my $expected = pop @line_elements;
		my $observed = pop @line_elements;
		my $description = join("\t", @line_elements);
		
		$expected = $observed / $expected;
		$vd_expected = $observed / $vd_expected;
		
		if( exists $interaction_counts{$description} ){
			${ $interaction_counts{$description} }[0] += $observed;
			${ $interaction_counts{$description} }[1] += $expected;
			${ $interaction_counts{$description} }[2] += $vd_expected;	
		}else{
			${ $interaction_counts{$description} }[0] = $observed;
			${ $interaction_counts{$description} }[1] = $expected;
			${ $interaction_counts{$description} }[2] = $vd_expected;
		}
	}
	close ANACONDA or die "Could not close filehandle on '$file' : $!";
}	

#Print out results
print OUT $header;
foreach my $interaction (keys %interaction_counts){
	my $total_observed = ${ $interaction_counts{$interaction} }[0];
	my $total_expected = ${ $interaction_counts{$interaction} }[1];
	my $total_vd_expected = ${ $interaction_counts{$interaction} }[2];
	
	my $total_obs_exp = $total_observed / $total_expected;
	my $total_obs_vd_exp = $total_observed / $total_vd_expected;
	
	print OUT "$interaction\t$total_observed\t$total_obs_exp\t$total_obs_vd_exp\n";	
}

close OUT or die "Could not close filehandle on '$outfile' : $!";

print "Combined files\n";

exit (0);
