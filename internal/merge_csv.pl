#!/bin/perl -w

use strict;
use warnings;

#Config
########

our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";

if( (scalar @ARGV != 3) || (! -f $ARGV[0]) || (! -f $ARGV[1]))
{
    die "Usage: script csv1 csv2 outputfile\n\n This script merges two csv files with headers by the first column. It assumes that the two have the same ids\n";
}

my $csv1=$ARGV[0];
my $csv2=$ARGV[1];
my $outputfile=$ARGV[2];

open(my $IFILE1, $csv1) or die "The file $csv1 cannot be read\n";
open(my $IFILE2, $csv2) or die "The file $csv1 cannot be read\n";

my @content1=<$IFILE1>;
my @content2=<$IFILE2>;

close $IFILE1;
close $IFILE2;

my $header1=splice(@content1,0,1);
my $header2=splice(@content2,0,1);
my $oheader;
chomp($header1);
chomp($header2);

my %parsed_csv1=%{parse_content(@content1)};
my %parsed_csv2=%{parse_content(@content2)};

if ((split($FS,$header1))[0] ne (split($FS,$header2))[0])
{
    die "The first column of the two csv files does not have the same name, make sure they are compatible and use the same name to use this script\n";
}

my @temp_header=split($FS,$header2);
shift(@temp_header);
$oheader=$header1.",".join($OFS,@temp_header);

open(my $OFILE, ">$outputfile") or die "The output file $outputfile cannot be open to be written\n";

print($OFILE "$oheader\n");

foreach my $key (keys %parsed_csv1)
{
    if (defined $parsed_csv2{$key})
    {
        #print("DEBUG $parsed_csv1{$key}\n");
        print($OFILE $key,$OFS,$parsed_csv1{$key},$OFS,$parsed_csv2{$key},"\n");
        delete($parsed_csv2{$key});
    }
    else
    {
        die "The id $key is not present in the second CSV. Make sure that the CSV files are compatible before using this script\n";
    }
}

if(scalar(keys %parsed_csv2) >=1)
{
    die ("The ids ".join($OFS,keys %parsed_csv2)." are not present in the first CSV. Make sure that the CSV files are compatible before using this script\n");
}
exit;

#Functions
###########

sub parse_content
{
    my %hash;
    my @temp;
    my @content=@_;
    my $key;
    foreach my $line (@content)
    {
        chomp($line);
        #print("DEBUG $line\n");
        @temp=split($FS,$line);
        $key=$temp[0];
        shift(@temp);
        $hash{$key}=join($OFS,@temp);
        #print("DEBUG $key, $hash{$key}\n");
    }
    return \%hash;
}
