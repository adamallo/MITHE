#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Cwd qw(cwd);

#Configuration variables
my $OFS="\t";
my $usage="$0 -d main_directory -o output_file [--covB covBprefix --paf PAFprefix]\n";

#Input variables
my $main_dir="";
my $covB="";
my $PAF="";
my $outfile="";
my $help;

##Getopt
######################################################
(! GetOptions(
    'main_directory|d=s' => \$main_dir,
    'output_file|o=s' => \$outfile,
    'covB|b=s' =>\$covB ,
    'PAF|p=s' =>\$PAF ,
    'help|h' => \$help,
                )) or (($outfile eq "") || ($main_dir eq "")  || $help) and die $usage;


open(my $OUTPUT,">$outfile") or die "Impossible to open the output file $ARGV[1]\n";

#my $output2name=$outfile;
#$output2name=~s/(\.[^.]+)$/_1col$1/g;
#open(my $OUTPUT2,">$output2name") or die "Impossible to open the output file $output2name\n";

#print($OUTPUT "Patient${OFS}Chr${OFS}Start${OFS}End${OFS}REF${OFS}ALT${OFS}Common${OFS}Common_pos${OFS}Common_type${OFS}Chr${OFS}Start${OFS}End${OFS}REF${OFS}ALT${OFS}A${OFS}A_pos${OFS}A_type${OFS}Chr${OFS}Start${OFS}End${OFS}REF${OFS}ALT${OFS}B${OFS}B_pos${OFS}B_type\n");

#print($OUTPUT2 "Patient${OFS}Shared${OFS}Chr${OFS}Start${OFS}End${OFS}REF${OFS}ALT${OFS}Name${OFS}Pos${OFS}Type\n");
print($OUTPUT "Patient${OFS}Shared${OFS}Chr${OFS}Start${OFS}End${OFS}REF${OFS}ALT${OFS}Name${OFS}Pos${OFS}Type\n");

chdir($main_dir);

opendir(my $DH,".");

my @dirs = grep {-d "\./$_" && ! /^\.{1,2}$/} readdir($DH);

closedir($DH);

foreach my $dir (@dirs)
{
    #print("DEBUG: $dir\n");
    my $name=basename($dir);

    my $filename="$dir/vcfdict.csv";
#    if (-f $filename)
#    {
#        print("DEBUG: opening ".cwd."/$filename\n");
#    }
#    else
#    { 
#        print("DEBUG: opening unexisting file ".cwd."/$filename\n");
#    }
    open(my $IFILE,$filename) or die "Error opening $filename\n";
    my @listfiles=<$IFILE>;
    close($IFILE);

    my @temp;
    my @common_files;
    my @da_files;
    my @db_files;

    for my $line (@listfiles)
    {
        chomp($line);
        #print("DEBUG: $line\n");
        if ($line =~ m/filt${covB}NABU${PAF}#.*common\.vcf/)
        {
            @temp=split(",",$line);
            push(@common_files,"$dir/".$temp[1].".annotated.variant_function");
           #print("\tDEBUG: pushing $temp[1] in common files\n");
        }
        elsif( $line =~ m/Afilt${covB}NAB${PAF}#.*different\.vcf/)
        {
            @temp=split(",",$line);
            push(@da_files,"$dir/".$temp[1].".annotated.variant_function");
            #print("\tDEBUG: pushing $temp[1] in da files\n");
        }
        elsif( $line =~ m/Bfilt${covB}NAB${PAF}#.*different\.vcf/)
        { 
            @temp=split(",",$line);
            push(@db_files,"$dir/".$temp[1].".annotated.variant_function");
            #print("\tDEBUG: pushing $temp[1] in db files\n");
        }
    }

#    my @common_files=glob("$dir/filt${covB}NABU#*common.vcf.annotated.variant_function");
    my $common_file=$common_files[0];
    my $common_file_aux=$common_file;
    $common_file_aux=~s/variant_function/exonic_variant_function/;

    #print("DEBUG: common files @common_files\n");
    if (scalar @common_files != 1)
    {
        die "More than one files detected with the format filt${covB}NABU${PAF}#*common.vcf.annotated.variant_function. This script is intended to analyse the output of only one filter\n";
    }

#   my @da_files=glob("$dir/Afilt${covB}NAB#*different.vcf.annotated.variant_function");
    my $da_file=$da_files[0];
    my $da_file_aux=$da_file;
    $da_file_aux=~s/variant_function/exonic_variant_function/;

    #print("DEBUG: privateA files @da_files\n");
    if (scalar @da_files != 1)
    {
        die "More than one files detected with the format Afilt${covB}NAB${PAF}#*different.vcf.annotated.variant_function. This script is intended to analyse the output of only one filter\n";
    }

#   my @db_files=glob("$dir/Bfilt${covB}NAB#*different.vcf.annotated.variant_function");
    my $db_file=$db_files[0];
    my $db_file_aux=$db_file;
    $db_file_aux=~s/variant_function/exonic_variant_function/;

    #print("DEBUG: privateB files @db_files\n");
    if (scalar @db_files != 1)
    {
        die "More than one files detected with the format Bfilt${covB}NAB${PAF}#*different.vcf.annotated.variant_function. This script is intended to analyse the output of only one filter\n";
    }

    open(my $COMMON,$common_file) or die "Impossible to open the input file $common_file\n";
    open(my $COMMON_AUX,$common_file_aux) or die "Impossible to open the input file $common_file_aux\n";
    open(my $DA, $da_file) or die "Impossible to open the input file $da_file\n";
    open(my $DA_AUX, $da_file_aux) or die "Impossible to open the input file $da_file_aux\n";
    open(my $DB, $db_file) or die "Impossible to open the input file $db_file\n";
    open(my $DB_AUX, $db_file_aux) or die "Impossible to open the input file $db_file_aux\n";

    my @c_content=<$COMMON>;
    my @c_content_aux=<$COMMON_AUX>;
    my @da_content=<$DA>;
    my @da_content_aux=<$DA_AUX>;
    my @db_content=<$DB>;
    my @db_content_aux=<$DB_AUX>;

    my @lengths=(scalar @c_content,scalar @da_content,scalar @db_content);
    @lengths=sort{$a <=> $b} @lengths;
   
    close($COMMON);
    close($DA);
    close($DB);
    close($COMMON_AUX);
    close($DA_AUX);
    close($DB_AUX);
    my @temp_data;
    my $type="";
    my @temp2=();

    for (my $i=0; $i<$lengths[2];++$i)
    {
#        my @out_line=("$name");
        if(scalar(@c_content) <= $i)
        {
            @temp_data=("","","","","","","","");
        }
        else
        {
            @temp_data=split("\t",$c_content[$i]);
            if($temp_data[0]=~m/^exonic$/)
            {
                foreach my $line (@c_content_aux)
                {
                    my $n=$i+1;
                    if($line=~/^line$n/)
                    {
                        @temp2=split("\t",$line);
                        $type=$temp2[1];
                        last;
                    }
                }
            }
            else
            {
                $type="";
            }
            $type=~s/$OFS/ /g;
            @temp_data=(@temp_data[2..6],$temp_data[1],$temp_data[0],$type);
            #print($OUTPUT2 print_array(($name,"Common",@temp_data)),"\n");
            print($OUTPUT print_array(($name,"Common",@temp_data)),"\n");
        }
        #print("DEBUG: @temp_data\n");
        #push(@out_line,@temp_data);
        #print("DEBUG: @out_line\n");
        if(scalar(@da_content) <= $i)
        {
            @temp_data=("","","","","","","","");
        }
        else
        {
            @temp_data=split("\t",$da_content[$i]);
           if($temp_data[0]=~m/^exonic$/)
            {
                foreach my $line (@da_content_aux)
                {
                    my $n=$i+1;
                    if($line=~/^line$n/)
                    {
                        @temp2=split("\t",$line);
                        $type=$temp2[1];
                        last;
                    }
                }
            }
            else
            {
                $type="";
            }

            $type=~s/$OFS/ /g;
            @temp_data=(@temp_data[2..6],$temp_data[1],$temp_data[0],$type);
            #print($OUTPUT2 print_array(($name,"A",@temp_data)),"\n");
            print($OUTPUT print_array(($name,"A",@temp_data)),"\n");
        }
        #push(@out_line,@temp_data);
        
        if(scalar(@db_content) <= $i)
        {
            @temp_data=("","","","","","","","");
        }
        else
        {
            @temp_data=split("\t",$db_content[$i]);
            if($temp_data[0]=~m/^exonic$/)
            {
                foreach my $line (@db_content_aux)
                {
                    my $n=$i+1;
                    if($line=~/^line$n/)
                    {
                        @temp2=split("\t",$line);
                        $type=$temp2[1];
                        last;
                    }
                }
            }
            else
            {
                $type="";
            }
            $type=~s/$OFS/ /g;
            @temp_data=(@temp_data[2..6],$temp_data[1],$temp_data[0],$type);
            #print($OUTPUT2 print_array(($name,"B",@temp_data)),"\n");
            print($OUTPUT print_array(($name,"B",@temp_data)),"\n");
        }
        #print("DEBUG: @temp_data\n");
        #push(@out_line,@temp_data);
#        print("DEBUG: @out_line\n");
        #print($OUTPUT print_array(@out_line),"\n");
    }
}

close($OUTPUT);
#close($OUTPUT2);
exit;

sub print_array
{
    my @array=@_;
    my $output="";
    foreach my $element (@array)
    {
        $output="$output$element$OFS";
    }
    chop($output);
    return($output);
}
