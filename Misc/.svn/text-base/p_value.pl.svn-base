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



print "Assigning P-values to data:\n";


foreach my $file (@files){
	print "\t$file\n";


	#Determine the number of times a feature is present in the dataset
	if($file =~ /\.gz$/){
		open (ANACONDA, "zcat $file |") or die "Could not read file '$file' : $!";
	}else{
		open(ANACONDA, '<', $file) or die "Could not open '$file' : $!";
	}
	scalar <ANACONDA>;    #Ignore headers
	
	my %feature_count;    # %{feature} = count	
	my $total_interactions = 0;
	while(<ANACONDA>){
		my $line = $_;
		chomp $line;
		my (undef, undef, undef, undef, $feature1, undef, undef, undef, undef, $feature2, $frequency) =  split(/\t/, $line);

		$feature_count{$feature1} += $frequency;
		$feature_count{$feature2} += $frequency;
		$total_interactions += $frequency;
	}

	close ANACONDA or die "Could not close filehandle on '$file' : $!";


	#Create a modified version of the file using with Poisson-derived P-value
	if($file =~ /\.gz$/){
		open (ANACONDA2, "zcat $file |") or die "Could not read file '$file' : $!";
	}else{
		open(ANACONDA2, '<', $file) or die "Could not open '$file' : $!";
	}

	my $outfile = "$file.pvalue.txt.gz";
	open(OUT, "| gzip -c - > $outfile") or die "Could not write to '$outfile' : $!";

	my $header = scalar <ANACONDA2>;
	chomp $header;
	$header .= "\tP-value(Poisson)\n";
	print OUT $header;

	while(<ANACONDA2>){
		my $line = $_;
		chomp $line;
		my (undef, undef, undef, undef, $feature1, undef, undef, undef, undef, $feature2, $frequency) =  split(/\t/, $line);

		my $feature1_count = $feature_count{$feature1}; 
		my $feature2_count = $feature_count{$feature2};

		my ($n1, $n3);   #n1 is larger than n3 (note called n3 here, not n2)
		if( $feature_count{$feature1} > $feature_count{$feature2} ){
			$n1 = $feature_count{$feature1};
			$n3 = $feature_count{$feature2};
		}else{
			$n1 = $feature_count{$feature2};
			$n3 = $feature_count{$feature1};
		}

		if( ($n1 / $total_interactions) > 0.05 ){   
			print OUT "$line\tNA\n";   #P-value could not be determined since n1 is too large 
		}else{
			my $lambda = ( ($n1 * $n3) / ( (2 * $total_interactions) - $n1 ) );
			my $p_value = poisson ($frequency, $lambda);
			print OUT "$line\t$p_value\n";
		}
	}

	close ANACONDA2 or die "Could not close filehandle on '$file' : $!";
	close OUT or die "Could not close filehandle on '$outfile' : $!";

}	


print "Assigned P-values\n";

exit (0);



##########################################################
#Subroutines
##########################################################

#Subroutine poisson:
#P(k \text{ events in interval}) = \frac{\lambda^k e^{-\lambda}}{k!}
#Takes k and lambda and input, returns p
sub poisson{
	my ($k, $lambda) = @_;
	$k = int($k);    #Should only be an integer anyway
	my $k_factorial = factorial($k);
	my $e = 2.7182818284;

	my $p = ( ($lambda ** $k) * ($e ** (-$lambda)) ) / $k_factorial;
	$p = printf("%.5f", $p);       # Rounded to 6 d.p.
	return $p;
}

#Subroutine factorial
#Takes an integer n and returns n!
sub factorial {
	my $n = $_[0];	
	my $f = 1;
	my $i = 1;

	$f *= ++$i while $i < $n;
	return $f;
}
















