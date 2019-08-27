#!/usr/bin/perl

#############################################
#Takes the data output from the complexes file
#and returns the data, filtering out reads
#which do not pass the barcode interactions
#threshold
#############################################

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;


#Option variables
my $help;
my $min;

#######################
#Get user-supplied parameters
my $config_result = GetOptions(
    "help" => \$help,   
    "minimum=i"    => \$min,
);

die "Could not parse options.\n" unless ($config_result);

if ( $help ) {
    print while (<DATA>);
    exit(0);
}

unless(@ARGV){
	die "Please specify a file to process.\n";
}

$min = 2 unless(defined $min);

print "Filtering out reads where N < $min\n";

foreach my $file (@ARGV){
	print "Processing '$file'\n";

	#Count reads
	my %counter;    #%{id} = count;
	my %summary;
	my @summary_categories = ( 'Reads_Processed', 'Reads_Kept', 'Reads_Discarded', 'Percentage_Reads_Kept' );
	foreach my $category (@summary_categories){
		$summary{$category} = 0;
	}
	
	if( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }

	while(<IN>){
		my $line = $_;
		chomp $line;		
		my ($id) = split(/\t/, $line);
		$counter{$id}++;
	}	
	close IN or die "Could not close '$file' : $!";
	
	if( $file =~ /\.gz$/ ) {
        open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
        open( IN, $file ) or die "Could not read '$file' : $!";
    }
	
	#Print out results
	my $outfile = "$file.filter_N" . "$min.txt";
	open(OUT, '>', $outfile) or die "Could not write to '$outfile' : $!";
	
	while(<IN>){
		my $line = $_;
		chomp $line;
		my ($id) = split(/\t/, $line);
		$summary{Reads_Processed}++;
		
		if($counter{$id} >= $min){
			print OUT "$line\n"; 
			$summary{Reads_Kept}++;
		}else{
			$summary{Reads_Discarded}++;
		}
	}	
	close IN or die "Could not close '$file' : $!";
	close OUT or die "Could not close '$outfile' : $!";
	
	#Print summary results
	$summary{Percentage_Reads_Kept} = calc_perc($summary{Reads_Kept}, $summary{Reads_Processed}, 2 );
	my $summary_file = "$file.group_N" . $min . '_summary.txt';
	open(SUMMARY, '>', $summary_file) or die "Could not open '$summary_file' : $!";
	print SUMMARY "File";
	foreach my $category (@summary_categories){
		print SUMMARY "\t$category";
	}
	print SUMMARY "\n";
	
	print SUMMARY "$file";	
	foreach my $category (@summary_categories){
		print SUMMARY "\t$summary{$category}";
	}
	print SUMMARY "\n";
	close SUMMARY or die "Could not close '$summary_file' : $!";	
}

print "Processing Complete.\n";

exit (0);



##############################################################################
#Subroutines
##############################################################################


######################
#Subroutine: calc_perc
#Receives a number and a total and returns the perentage value
#Optional argument: decimal places of the output
#Subroutine rounds following the sprintf rounding protocol 
sub calc_perc {
    my ( $n, $total, $dp ) = @_;
    
    if(defined $dp){
        $dp = abs( int($dp) );
    }else{
        $dp = 2;
    }
    
    return 'NA' unless(defined $n and defined $total);   #Avoid initialisation error
    return 'NA' if($total == 0);    #Avoid division by zero error
    
    my $pc = 100 * $n / $total;
    my $pc_string = '%.' . $dp . 'f';

    $pc = sprintf("$pc_string", $pc);

    
    return $pc;
}








__DATA__

Takes the data output from the complexes pipeline  and returns the data, filtering out reads
which do not pass the barcode group size 
threshold

OPTIONS

--help         Print program help and exit
--minimum      Minimum size of complex (default=2)

Steven Wingett, Babraham Institute, Cambridge, UK (steven.wingett@babraham.ac.uk)
