#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd;
use File::Basename;

##Configuration variables
######################################################

our $sep_param="#";
our $sep_value=":";
our $OFS=",";
our $FS=",";
#our $n_cores=1;

######################################################

##IO Variables
######################################################
my $input_file="";
my $output_file="";
my ($remove_all_but,$remNAB,$remcovB,$rempopAF)=("",0,0,0); ##By default we keep all of them
#my $original_dir="";

#Flags
my $help;
my $usage="Usage: $0 [options] -i input_file -o output_file \n\nOptions:\n--------\n\t--remove_all_but: List of comma-separated column names to keep\n\t--remove_NAB_collapse: Remove the NAB section and collapse the resulting equivalent filtering conditions\n\t--remove_covB_collapse: Remove the covB section and collapse the resulting equivalent filtering conditions\n\t--remove_popAF_collapse: Remove the popAF section and collapse the resulting equivalent filtering conditions\n\n";
######################################################

######################################################
##MAIN
######################################################

##Getopt
######################################################
((! GetOptions(
    'remove_all_but=s' => \$remove_all_but,
	'input_file|i=s' => \$input_file,
	'output_file|o=s' => \$output_file,
	'remove_NAB_collapse=i' => \$remNAB,
    'remove_covB_collapse=i' => \$remcovB,
    'remove_popAF_collapse=i' => \$rempopAF,
    #'n_cores=i' => \$n_cores,
    'help|h' => \$help,
                )) or (($output_file eq "") || ($input_file eq "") || $help)) and die $usage;

##IO
open(my $IFILE, $input_file);
open(my $OFILE, ">$output_file");

my @content= <$IFILE>;

close($IFILE);

##Get Conditions header
my @temp;
my @conditions_header=@{get_condition_parameters($content[1])};
##Add to the header the needed columns and store their indexes
my @indexes;
my @original_header=split($FS,$content[0]);
chomp(@original_header);
my $n_columns=scalar @original_header;

if($remove_all_but eq "")
{
	@indexes=1 .. $n_columns-1;
}
else
{
	my %avaliable_columns;
	@avaliable_columns{@original_header}=(0 .. scalar @original_header); ##This initiates avaliable_columns with temp as keys and its index as values
	my @columns=split($FS,$remove_all_but);
	chomp(@columns);
	for (my $i=0; $i<scalar(@columns); ++$i)
	{
		if(exists $avaliable_columns{$columns[$i]})
		{
			push(@indexes,$avaliable_columns{$columns[$i]});
			delete($avaliable_columns{$columns[$i]});
		}
		else
		{
			die "The column $columns[$i] is not present in the original file. Check your input options\n";
		}
	}
}

#Print the new header
print($OFILE join($OFS,@conditions_header,@original_header[@indexes]),"\n"); ##Prints the header

#print("DEBUG: ", join($OFS,@header,@original_header[@indexes]), "\n");

##Row loop
####################
my %out_hash; ##Using a hash equivalent conditions are automatically collapsed

for (my $row=1;$row<scalar @content; ++$row)
{
	chomp($content[$row]);
    my @columns=split($FS,$content[$row]);
	my $condition=join($OFS,@{get_condition_values($columns[0])});
	$out_hash{$condition}=join($OFS,@columns[@indexes]);
	##print($OFILE join($OFS,$condition,@columns[@indexes]),"\n");
}

foreach my $key (keys %out_hash)
{
	print($OFILE join($OFS,$key,$out_hash{$key}),"\n");
}

close($OFILE);

sub get_condition_parameters
{
	my ($textheader)=@_;
	my @header;
	$textheader=(split($FS,$textheader))[0];
	$textheader=~s/--//g;##Remove the -- to make them more human-readable
	$textheader=~s/$sep_value/ /g;
	$textheader=~s/$sep_param/ /g;
	$textheader=~tr/=/ /;
	@temp=split(" ",$textheader);
	#@temp=split($sep_value,join($sep_value,split("=",join($sep_value,split($sep_param,(split($FS,$content[1]))[0])))));
	my $concat="";
	my $wait=0;

	for (my $i=0; $i<scalar @temp; $i+=2)
	{ 
		if($temp[$i] eq "covB")
        {
            --$i;
            $concat="covB";
            if($remcovB==1)
            {
                $wait=1;
            }
            next;
        }
		elsif($temp[$i] eq "NAB")
		{
            $wait=0;
			--$i;
			$concat="NAB";
			if($remNAB==1)
			{
                $wait=1;
			}
			next;
		}
        elsif($temp[$i] eq "PAF")
        {
            $wait=0;
            --$i;
            $concat="PAF";
            if($rempopAF==1)
            {
				last; ##If I added more of these special conditions (like covB or NAB) I would just need to put them in order modifying the wait. The last can exit earlier with last, like here
            }
            next;
        }
        elsif($wait == 1)
        {
            next;
        }
		else
		{
			push(@header,"$concat$temp[$i]");
		}
	}
	
	return \@header;
}
sub get_condition_values
{
	my ($textheader)=@_;
	my @header;
	$textheader=~s/$sep_value/ /g;
	$textheader=~s/$sep_param/ /g;
	$textheader=~tr/=/ /;
	@temp=split(" ",$textheader);
	#@temp=split($sep_value,join($sep_value,split("=",join($sep_value,split($sep_param,(split($FS,$content[1]))[0])))));
	##my $concat="";
	my $wait=0;

	for (my $i=1; $i<scalar @temp; $i+=2)
	{ 
		if($temp[$i-1] eq "covB")
        {
            --$i;
            #$concat="covB";
            if($remcovB==1)
            {
                $wait=1;
            }
            next;
        }
		elsif($temp[$i-1] eq "NAB")
		{
            $wait=0;
			--$i;
            #$concat="NAB";
			if($remNAB==1)
			{
                $wait=1;
			}
			next;
		}
        elsif($temp[$i-1] eq "PAF")
        {
            $wait=0;
            --$i;
            #$concat="PAF";
            if($rempopAF==1)
            {
				last; ##If I added more of these special conditions (like covB or NAB) I would just need to put them in order modifying the wait. The last can exit earlier with last, like here
            }
            next;
        }
        elsif($wait == 1)
        {
            next;
        }
		else
		{
            #push(@header,"$concat$temp[$i]");
			push(@header, $temp[$i]);
		}
	}
	return \@header;
}


