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
#Perl script takes a tab-delimited data file of
#Read_ID	
#Barcode_ID	
#Chromosome	
#Feature_Start	
#Feature_End	
#Feature_Strand	
#Feature_ID	
#Feature_Name
#
#Returns a file of virtual di-tags
#Chromosome	
#Feature_Start	
#Feature_End	
#Feature_Strand	
#
##############################################################################

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/../";
use module_anaconda;

use Data::Dumper;

#Check input ok
die "Please specify files to process.\n" unless (@ARGV);
die "Please adjust configuration.\n" unless( check_files_exist(\@ARGV, 'EXISTS') );

foreach my $file (@ARGV){
    
    print "Processing $file\n";

	if( $file =~ /\.gz$/ ) {
    	open( IN, "zcat $file |" ) or die "Couldn't read '$file' : $!";
	} else {
    	open( IN, '<', $file ) or die "Could not read '$file' : $!";
	}
    
    my $ditags_outfile = $file;
    $ditags_outfile =~ s/\.simulation_formatted\.txt\.gz$//;
    $ditags_outfile .= '.seqmonk_virtual_ditags.txt.gz';

   	open(DITAGS, "| gzip -c - >  $ditags_outfile") or die "Could not write to '$ditags_outfile' : $!";
    open (CIS_TRANS, '>', "$file.cis_trans_by_group.txt") or die "Could not read '$file.cis_trans_by_group.txt' : $!";
    print CIS_TRANS "Group_Size\tPercentage_Trans\n";
    
    my %groupedReads;    #%{group} = @(reads)	
    scalar <IN>;    #Ignore header
    while(<IN>){
		my $line = $_;
		chomp $line;
		my (undef, $group, $csome, $start, $end, $strand) = split(/\t/, $line);
		push ( @{ $groupedReads{$group} }, "$csome\t$start\t$end\t$strand");
    }
    
	#print Dumper \%groupedReads;
    
    #Create the virtual HiC dataset
    foreach my $group (keys %groupedReads){
		next unless( scalar(@{ $groupedReads{$group} }) > 1);    #Only process if Barcode group size > 1
		my @ditags = createDitags( @{ $groupedReads{$group} } );
		print DITAGS join("\n", @ditags);
		print DITAGS "\n";
	
		#Calculate the percentage trans for each group size
		my $group_size = scalar @{ $groupedReads{$group} };
		my $trans = calcTrans(@ditags);
		print CIS_TRANS "$group_size\t$trans\n";
    }
    
    close IN or die "Could not close '$file' : $!";
    close DITAGS or die "Could not close '$ditags_outfile' : $!";
    close CIS_TRANS or die "Could not close '$file.cis_trans_by_group.txt' : $!";
	
}

print "Processing complete.\n";

exit (0);





#######################################################################################
#Subroutines                                                                          #
#######################################################################################

#calcTrans
#Takes an array of di-tags, with paired reads on adjacent lines,
#and reports the percentage trans di-tags
sub calcTrans{
	my @ditags = @_;
	my $cis = 0;
	my $trans = 0;
	
	while(scalar @ditags > 0){
		my $read1 = shift @ditags;
		my $read2 = shift @ditags;
		
		my ($csome1) = split(/\t/, $read1);
		my ($csome2) = split(/\t/, $read2);
		
		$csome1 eq $csome2 ? $cis++ : $trans++;
	}
	
	my $percTrans = 100 * $trans / ($cis + $trans);
	return $percTrans;
}



#createDitags
#Takes an array of reads and creates an array of the read-pair combinations
#(reports read pairs once i.e. read1-read2 but not read2-read1 as well)
sub createDitags{
	my @reads = @_;
	my @ditags;
 
	while(scalar @reads > 1){
		for (my $i=1; $i < scalar @reads; $i++){
			push(@ditags, $reads[0]);
			push(@ditags, $reads[$i]);
		}	
		shift @reads;
	}
	#print Dumper \@ditags;
	return @ditags;
}
