#!/bin/perl -w

use warnings;
use strict;

#Configuration parameters
#########################

our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";

###CHECK ARGV
if(scalar @ARGV != 2)
{
    die "Usage: script inputfile outputfile\n";
}
my $inputfile=$ARGV[0];
my $outputfile=$ARGV[1];

open(my $IFILE, $inputfile) or die "Error reading the input file $inputfile\n";
my @content=<$IFILE>;
close($IFILE);


#Hashes, keys=filtering_condition, value=tstv
my %Aconds;
my %Bconds;
my %filtUconds;
my %filtIconds;
my %filtNUconds;
my %filtNIconds;
my %filtcovBNUconds;
my %filtcovBNIconds;
my %fullUconds;
my %fullIconds;
my %fullcovBUconds;
my %fullcovBIconds;
my %fullcovBpopAFUconds;
my %fullcovBpopAFIconds;

my @temp;

foreach my $line (@content)
{
    chomp($line);
    @temp=split($FS,$line);
    #print("DEBUG: $temp[0], $temp[1]\n");
    if($temp[0]=~m/different/)
    {
        next;
    }
    if($temp[0]=~s/^filtNABU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $fullUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtNABI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $fullIconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtcovBNABUPAF${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $fullcovBpopAFUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtcovBNABIPAF${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $fullcovBpopAFIconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtcovBNABU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $fullcovBUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtcovBNABI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $fullcovBIconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^A${sep_param}//)
    {
        $Aconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^B${sep_param}//)
    {
        $Bconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtNU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtNUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtNI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtNIconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtNcovBU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtcovBNUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtNcovBI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtcovBNIconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtU${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtUconds{$temp[0]}=$temp[1];
    }
    elsif($temp[0]=~s/^filtI${sep_param}//)
    {
        $temp[0]=~s/_common//;
        $filtIconds{$temp[0]}=$temp[1];
    }

}

open(my $OFILE, ">$outputfile") or die "The file $outputfile can't be opened, check your input options\n";
print($OFILE join($OFS,("Condition","TsTv_A","TsTv_B","TsTv_filtU","TsTv_filtI","TsTv_filtNU","TsTv_filtNI","TsTv_filtNcovBU","TsTv_filtNcovBI","TsTv_filtNABU","TsTv_filtNABI","TsTv_filtcovBNABU","TsTv_filtcovBNABI","TsTv_filtcovBNABUPAF","TsTv_filtcovBNABIPAF")),"\n");

###Pending: I should make sure that I and U have the same keys for every filtering level

my $filtcond;
my $covBcond;
my $NABcond;
my $popAFcond;

foreach my $fullcond (keys %fullcovBpopAFUconds)
{
    ($filtcond,$covBcond,$NABcond,$popAFcond)=($fullcond=~m/^(.*?)(${sep_param}covB${sep_param}.*?)(${sep_param}NAB${sep_param}.*?)(${sep_param}PAF${sep_param}.*?)$/);
    print($OFILE join($OFS,($fullcond,$Aconds{$filtcond},$Bconds{$filtcond},$filtUconds{$filtcond},$filtIconds{$filtcond},$filtNUconds{$filtcond},$filtNIconds{$filtcond},$filtcovBNUconds{"$filtcond$covBcond"},$filtcovBNIconds{"$filtcond$covBcond"},$fullUconds{"$filtcond$NABcond"},$fullIconds{"$filtcond$NABcond"},$fullcovBUconds{"$filtcond$covBcond$NABcond"},$fullcovBIconds{"$filtcond$covBcond$NABcond"},$fullcovBpopAFUconds{$fullcond},$fullcovBpopAFIconds{$fullcond})),"\n");

}
close($OFILE);
exit;
