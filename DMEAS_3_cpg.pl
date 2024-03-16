#!/usr/local/bin/perl
use POSIX;
use strict;
##########################3

######                           ###
  ##     #    #                     #   
  ##    ###   #      ###           #
  ##     #          #           ###
  ##     #           ###           #     and 5 reads
######   ##          #         ####
##################################################################################
# Description: 
# This perl script is devoleped as a complementary to DMEAS to calculate methylation entropy in Linux/Unix workstation. 
# This script can handle both genome-wide methylation data and locus-specific methylation data as shown in DMEAS.
# This script supports in both batch-entry mode (with specifying the input files folder) and single-entry mode (with specifying the input file). 
##################################################################################
#
# As DMEAS supports both multiple-entry mode and single-entry mode, specific commands need to be executed. For multiple-entry mode, the corresponding folder/directory path (which contain the methylation data) needs to be specified, while for the single-entry mode, the corresponding file path needs to be specified.
# The typical command-line could look like as follows:
# perl DMEAS.pl [options] input(s) output
# In the command-line, the input(s) and output must be provided whereas the argument [options] is optional to specify the usage of the script. 
# 
# Usage:
# perl DMEAS.pl [-d] [gw|ls] [-s|-f] [ms/ss|sf] input1 [input2] -o output_folder 
#
# The parameters of this perl script are as follows:
# 
# -d                 An parameter (combined with [gw|ls]) specifies the DNA methylation data type such as genome-wide methylation data and locus-specific methylation data. This parameter is optional.
# gw|ls              This optional parameter is only available when the parameter '-d' is used. Here 'gw' represents genome-wide methylation data and 'ls' represents locus-specific methylation data.
# -s|-f              An parameter (combined with [ms/ss|sf]) specifies the mode of input files such as batch-entry and single-entry. The parameter '-s' should be used if the genome-wide methylation data is input while '-f' should be used if the locus-specific methylation data is input. This parameter is optional.
# ms/ss|sf           This optional parameter is only available when the parameter '-s|-f' is used. When '-s' is used, the 'ms/ss' will be followed and will specify the mode of input files of genome-wide methylation data, i.e., 'ms' specifies the batch-entry mode while 'ss' specifies the single-entry mode. When '-f' is used, the 'sf' will be followed and will specify the mode of input files of locus-specific methylation data. Currently, only single-file input is supported by DMEAS for locus-specific methylation data.
# input1             This argument should be provided for both genome-wide and locus-specific methylation data. Regarding to genome-wide methylation data (as '-d gw'), there are two options as indicated by '-s ms/ss'. That is, the full path to the specific input file containing the genome-wide methylation data should be provided as single-entry mode is used ('-s ss') while the full path to the folder containing the multiple entries of genome-wide methylation data should be provided as batch-entry mode is used ('-s ms'). Regarding to locus-specific methylation data (as '-d ls'), only the full path to the file containing locus-specific methylation data should be provided (as '-f sf'). 
# input2             This input file is only available when the parameter '-d' is used. As '-d' will specify the genome-wide methylation data, the CpG genome position information (which is saved in input2) should be provided.
# -o                 A required parameter indicating the output.
# output_folder      The full path to the folder which will contain all the output results from DMEAS.
# 
# Note: the order of parameters should be strictly followed as provided above.
#
# For the genome-wide methylation data with multiple files (batch-entry), the command-line looks like as follows:
#   perl DMEAS.pl -d gw -s [ss|ms] folder_path_bis folder_path_gp -o output_path
#
# For example:
#   perl DMEAS.pl -d gw -s ms ./test/bismark/ ./test/genome_position/ -o ./output
#   perl DMEAS.pl -d gw -s ss ./test/bismark/21.bis ./test/genome_position/21.gp -o ./output
# Here, the argument folder_path_bis represents a folder path contains files with ".bis" (no other files allowed in the given folder); the argument folder_path_gp represents a folder path contain the CpG genome position files with ".gp" (no other files as well). The ".gp" must correspond to the files with ".bis" accordingly. The argument output_path represents the output's folder path where the user would like to save all the output results. 
#
# The default parameter of sample types for genome-wide methylation data ("-d gw") is "-s ms", thus it could be omitted unless single sample bismark methylation data is specified. The command-line looks like as follows:
#   (1) perl DMEAS.pl folder_path_bis folder_path_gp -o output_path
#   (2) perl DMEAS.pl -s [ss|ms] file_path_bis file_path_gp -o output_path
#
# For example:
#   perl DMEAS.pl ./test/bismark/ ./test/genome_position/ -o ./output
#   perl DMEAS.pl -s ss ./test/bismark/21.bis ./test/genome_position/21.gp -o ./output
#   perl DMEAS.pl -s ms ./test/bismark/ ./test/genome_position/ -o ./output
# In the first case (1), the DMEAS can automatically perform analysis on the genome-wide methylation data with multiple chromosomes. In the second case (2), the argument "-s ss" specifies the single methylation data, thus the detailed files of file_path_bis (".bis") and file_path_gp (".gp") are required to be provided.
#
# 
# For the locus-specific methylation data, the usage could look like this:
#   perl DMEAS.pl -d ls -f sf file_path_txt -o output_path
#
# For example:
#   perl DMEAS.pl -d ls -f sf ./test/locus/locus_1.txt -o ./output
# Here, the argument file_path_txt represents a file path contains the locus-specific methylation files with ".txt". The command-line looks like as follows: 
#   (1) perl DMEAS.pl file_path_txt -o output_path
#   (2) perl DMEAS.pl -f sf file_path_txt -o output_path
# For example:
#   perl DMEAS.pl ./test/locus/locus_1.txt -o ./output
#   perl DMEAS.pl -f sf ./test/locus/locus_1.txt -o ./output
# Similarly, in the first case (1), DMEAS can automatically perform analysis on the locus-specific methylation data with single file by identifying the number of input command-line parameters. In the second case (2), the argument "-f" specifies the single methylation data for locus-specific methylation data analysis.
#
##################################################################################
# Error List:
#  error 01: Input arguments are too less or too many;
#  error 02: Input an invalid file path;
#  error 03: Input an invalid folder path;
#  error 04: Input arguments are invalid according to the rule of setting options;
#  error 05: The parameters' abbreviation is error;
#  error 06: The files with .gp saved genome position data isn't corresponding to the files with .bis saved methylation data from Bismark.
##################################################################################
# Warning List:
#  warning 01: Using the default parameters -d(ls) or -f(sf) for locus-specific methylation data;
#  warning 02: Using the default parameters -d(gw) or -s(ms) for genome-wide methylation data;
##################################################################################
# Version:1.01v.
##################################################################################

my $datatype="";
my $sampletype="";
my $filetype="";
my $bismark="";
my $genomewide="";
my $locusspecific="";
my $outpath="";

if(scalar(@ARGV)<3)
{
	print "Error Number 01.\n";
	print "Input arguments are too less.\n";
	print "When input three parameters, the default data type is Locus-Specific and the default type of input methylation data is single-entry.\n";
	print "When input four parameters, the default data type is Genome-wide and the default sample type is multi-entry.\n";
	print "Please refer user manual for more information.\n";
	exit(0);
}
elsif(scalar(@ARGV)==3)
{
	open(LOCUSSPECIF,"<$ARGV[0]") or die"Error Number 02.\nInput an invalid file path $ARGV[0] saved Locus-Specific methylation data.\n";
	close(LOCUSSPECIF);
  if($ARGV[1]!~/-o/i)
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
	if(!-d $ARGV[2])
	{
		mkdir $ARGV[2];
	}
	$outpath=$ARGV[2];
	print "Warning 01.\n";
	print "NOTE: here the default data type is Locus-Specific and the default type of saving methylation data is single file.\n";
	print "Please refer user manual for more information.\n";
	$datatype="ls";
	$filetype="sf";
	$locusspecific=$ARGV[0];
}
elsif(scalar(@ARGV)==4)
{
	opendir(BISMARK,"$ARGV[0]") or die"Error Number 03.\nInput an invalid folder path $ARGV[0] saved methylation data from Bismark.\n";
	closedir(BISMARK);
	opendir(GENOMEPOSITION,"$ARGV[1]") or die"Error Number 03.\nInput an invalid folder path $ARGV[1] saved genomic positions.\n";
	closedir(GENOMEPOSITION);
	if($ARGV[2]!~/-o/i)
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
	if(!-d $ARGV[3])
	{
		mkdir $ARGV[3];
	}
	$outpath=$ARGV[3];
	print "Warning 02.\n";
	print "NOTE: here the default data type is Genome-Wide and the default sample type is Multi-Sample.\n";
	print "Moreover, the two input parameters represent the folder saved methylation data from Bismark and the corresponding folder saved genomic positions respectively.\n";
	print "Please refer user manual for more information.\n";
	$datatype="gw";
	$sampletype="ms";
  $bismark=$ARGV[0];
  $genomewide=$ARGV[1];
}
elsif(scalar(@ARGV)==5)
{
	if($ARGV[0]=~/-f/i)
  {
  	#if($ARGV[1]=~/[m?]f/i)
		#{
  	#  opendir(LOCUSSPECIF,"$ARGV[2]") or die"Error Number 03.\nInput an invalid folder path $ARGV[2] saved Locus-Specific methylation data.\n";
	  #  closedir(LOCUSSPECIF);
    #}
    if($ARGV[1]=~/[s?]f/i)
    {
    	open(LOCUSSPECIF,"<$ARGV[2]") or die"Error Number 02.\nInput an invalid file path $ARGV[2] saved Locus-Specific methylation data.\n";
	    close(LOCUSSPECIF);
    }
    else
    {
      print "Error Number 04.\n";
    	print "Input arguments are invalid. \n";
      print "Please refer user manual for more information.\n";
      exit(0);
    }
    $locusspecific=$ARGV[2];
  }
  elsif($ARGV[0]=~/-d/i)
  {
    open(LOCUSSPECIF,"<$ARGV[2]") or die"Error Number 02.\nInput an invalid file path $ARGV[2] saved Locus-Specific methylation data.\n";
	  close(LOCUSSPECIF);
	  $locusspecific=$ARGV[2];
  }
  else
  {
    print "Error Number 04.\n";
    print "Input arguments are invalid. \n";
    print "Please refer user manual for more information.\n";
  	exit(0);
  }
  if($ARGV[3]!~/-o/i)
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
  if(!-d $ARGV[4])
	{
		mkdir $ARGV[4];
	}
	$outpath=$ARGV[4];
	if($ARGV[0]=~/-d/i)
	{
		print "Warning 01.\n";
	  print "NOTE: here the default type of saving methylation data is single file.\n";
	  print "Please refer user manual for more information.\n";
	  $datatype=$ARGV[1];
	  $filetype="sf";
	}
	elsif($ARGV[0]=~/-f/i)
	{
		print "Warning 01.\n";
		print "NOTE: here the default data type is Locus-Specific.\n";
		print "Please refer user manual for more information.\n";
		$datatype="ls";
	  $filetype=$ARGV[1];
	}
	else
	{
		print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
	}
}
elsif(scalar(@ARGV)==6)
{
	if($ARGV[0]=~/-s/i)
	{
  	if($ARGV[1]=~/[m?]s/i)
		{
  	  opendir(BISMARK,"$ARGV[2]") or die"Error Number 03.\nInput an invalid folder path $ARGV[2] saved methylation data from Bismark.\n";
  	  closedir(BISMARK);
  	  opendir(GENOMEPOSITION,"$ARGV[3]") or die"Error Number 03.\nInput an invalid folder path $ARGV[3] saved genomic positions.\n";
  	  closedir(GENOMEPOSITION);
    }
    elsif($ARGV[1]=~/[s?]s/i)
    {
      open(BISMARK,"<$ARGV[2]") or die"Error Number 02.\nInput an invalid file path $ARGV[2] saved methylation data from Bismark.\n";
  	  close(BISMARK);
  	  open(GENOMEPOSITION,"<$ARGV[3]") or die"Error Number 02.\nInput an file path $ARGV[3] saved genomic positions.\n";
  	  close(GENOMEPOSITION);
    }
    else
    {
    	print "Error Number 04.\n";
  		print "Input arguments are invalid. \n";
  		print "Please refer user manual for more information.\n";
  		exit(0);
    }
  	$bismark=$ARGV[2];
    $genomewide=$ARGV[3];
  }
  elsif($ARGV[0]=~/-d/i)
  {
  	opendir(BISMARK,"$ARGV[2]") or die"Error Number 03.\nInput an invalid folder path $ARGV[2] saved methylation data from Bismark.\n";
  	closedir(BISMARK);
  	opendir(GENOMEPOSITION,"$ARGV[3]") or die"Error Number 03.\nInput an invalid folder path $ARGV[3] saved genomic positions.\n";
  	closedir(GENOMEPOSITION);
  	$bismark=$ARGV[2];
    $genomewide=$ARGV[3];
  }
  else
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
  	exit(0);
  }
  if($ARGV[4]!~/-o/i)
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
  if(!-d $ARGV[5])
	{
		mkdir $ARGV[5];
	}
	$outpath=$ARGV[5];
	if($ARGV[0]=~/-d/i)
	{
		print "Warning 02.\n";
		print "NOTE: here the default sample type is Multi-Sample.\n";
		print "Please refer user manual for more information.\n";
		$datatype=$ARGV[1];
	  $sampletype="ms";
	}
	elsif($ARGV[0]=~/-s/i)
	{
		print "Warning 02.\n";
		print "NOTE: here the default data type is Genome-Wide.\n";
		print "Please refer user manual for more information.\n";
		$datatype="gw";
	  $sampletype=$ARGV[1];
	}
	else
	{
		print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
	}
}
elsif(scalar(@ARGV)==7)
{
	if($ARGV[1]=~/ls/i)
	{
	  #if($ARGV[3]=~/[m?]f/i)
		#{
  	#  opendir(LOCUSSPECIF,"$ARGV[4]") or die"Error Number 03.\nInput an invalid foler path $ARGV[4] saved Locus-Specific methylation data.\n";
	  #  closedir(LOCUSSPECIF);
    #}
    if($ARGV[3]=~/[s?]f/i)
    {
    	open(LOCUSSPECIF,"<$ARGV[4]") or die"Error Number 02.\nInput an invalid file path $ARGV[4] saved Locus-Specific methylation data.\n";
	    close(LOCUSSPECIF);
    }
    else
    {
    	print "Error Number 04.\n";
  		print "Input arguments are invalid. \n";
  		print "Please refer user manual for more information.\n";
  		exit(0);
    }
	  $locusspecific=$ARGV[4];
  }
  else
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
  if($ARGV[5]!~/-o/i)
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
  if(!-d $ARGV[6])
	{
		mkdir $ARGV[6];
	}
	$outpath=$ARGV[6];
	if($ARGV[0]=~/-d/i)
	{
		$datatype=$ARGV[1];
	}
	else
	{
		print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
	}
	if($ARGV[2]=~/-f/i)
	{
		$filetype=$ARGV[3];
	}
	else
	{
		print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
	}
}
elsif(scalar(@ARGV)==8)
{
	if($ARGV[1]=~/gw/i)
	{
		if($ARGV[3]=~/[m?]s/i)
		{
  	  opendir(BISMARK,"$ARGV[4]") or die"Error Number 03.\nInput an invalid folder path $ARGV[4] saved methylation data from Bismark.\n";
  	  closedir(BISMARK);
  	  opendir(GENOMEPOSITION,"$ARGV[5]") or die"Error Number 03.\nInput an invalid folder path $ARGV[4]saved genomic positions.\n";
  	  closedir(GENOMEPOSITION);
    }
    elsif($ARGV[3]=~/[s?]s/i)
    {
    	open(BISMARK,"<$ARGV[4]") or die"Error Number 02.\nInput an invalid file path $ARGV[4] saved methylation data from Bismark.\n";
  	  close(BISMARK);
  	  open(GENOMEPOSITION,"<$ARGV[5]") or die"Error Number 02.\nInput an invalid file path $ARGV[5] saved genomic positions.\n";
  	  close(GENOMEPOSITION);
    }
    else
    {
    	print "Error Number 04.\n";
  		print "Input arguments are invalid. \n";
  		print "Please refer user manual for more information.\n";
  		exit(0);
    }
	  $bismark=$ARGV[4];
    $genomewide=$ARGV[5];
  }
  else
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
  if($ARGV[6]!~/-o/i)
  {
  	print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
  }
  if(!-d $ARGV[7])
	{
		mkdir $ARGV[7];		
	}
	$outpath=$ARGV[7];
  if($ARGV[0]=~/-d/i)
	{
		$datatype=$ARGV[1];
	}
	else
	{
		print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
	}
	if($ARGV[2]=~/-s/i)
	{
	  $sampletype=$ARGV[3];
	}
	else
	{
		print "Error Number 04.\n";
		print "Input arguments are invalid. \n";
		print "Please refer user manual for more information.\n";
		exit(0);
	}
}
elsif(scalar(@ARGV)>=9)
{
	print "Error Number 01.\n";
  print "Input arguments are invalid. \n";
  print "Please refer user manual for more information.\n";
  exit(0);
}

if(!-d $outpath)
{
	if(scalar(split(/\\/,$outpath))>1 && scalar(split(/\//,$outpath))<=1)
	{
		my @tempArray=split(/\\/,$outpath);
		my $str=$tempArray[0];
		foreach my $i (1..scalar(@tempArray)-1)
		{
			$str.="/".$tempArray[$i];
			if(!-d $str)
			{
				mkdir $str;
			}
		}
	}
	elsif(scalar(split(/\\/,$outpath))<=1 && scalar(split(/\//,$outpath))>1)
	{
		my @tempArray=split(/\//,$outpath);
		if(scalar(@tempArray)!=0)
	  {
		  my $str=$tempArray[0];
		  foreach my $i (1..scalar(@tempArray)-1)
	  	{
		  	$str.="/".$tempArray[$i];
		  	if(!-d $str)
	  		{
		  		mkdir $str;
			  }
		  }
	  }
	  else
	  {
	  	print "error 03:\ninput an invalid folder path and please check the path $outpath.\n";
	  	exit(0);
	  }
	}
	else
	{
		print "error 03:\ninput an invalid folder path and please check the path $outpath.\n";
		exit(0);
	}
}
	
if($datatype=~/gw/i)
{
	if($sampletype=~/ms/i)
	{
		my %bismarkdata=();
		my %genomepositiondata=();
		my @tempAray=();

    opendir(DIR,"$bismark") or die"Error Number 03.\nInput an invalid folder saved methylation data from Bismark.\n";
    while(my $temp=readdir(DIR))
    {
    	if($temp=~/.+(\.bis)/)
    	{
    		@tempAray=split(/\./,$temp);
  	  	if(!$bismarkdata{$tempAray[0]})
   	  	{
	  	  	$bismarkdata{$tempAray[0]}=$bismark."/".$temp;
      	}
    	}
    }
    closedir(DIR);
    
    @tempAray=();
    opendir(DIR,"$genomewide") or die"Error Number 03.\nInput an invalid folder saved genomic positions.\n";
    while(my $temp=readdir(DIR))
    {
    	if($temp=~/.+(\.gp)/)
    	{
    		@tempAray=split(/\./,$temp);
  	  	if(!$genomepositiondata{$tempAray[0]})
   	  	{
	  	  	$genomepositiondata{$tempAray[0]}=$genomewide."/".$temp;
      	}
    	}
    }
    closedir(DIR);
    
    @tempAray=();
    foreach my $key (keys %bismarkdata)
    {
    	if(!$genomepositiondata{$key})
    	{
    		print "Error Number 06.\n";
        print "The files with .gp saved genome position data must be corresponding to the files with .bis saved methylation data from Bismark. \n";
        exit(0);
    	}
    	else
    	{
    		my $OUT=$outpath."/".$key;
    		&calmethynentropy($bismarkdata{$key},$genomepositiondata{$key},$OUT);
    	}
    }
	}
	elsif($sampletype=~/ss/i)
	{
		&calmethynentropy($bismark,$genomewide,$outpath);
	}
	else
	{
		print "Error Number 05.\n";
    print "Unexpected Bug. \n";
    exit(0);
	}
}
elsif($datatype=~/ls/i)
{
	#if($filetype=~/mf/i)
	#{
	#	my %locusspecificdata=();
	#	my @tempAray=();
	#	
	#	opendir(DIR,"$locusspecific") or die"Error Number 03.\nInput an invalid folder saved methylation data from Bismark.\n";
  #  while(my $temp=readdir(DIR))
  #  {
  #  	if($temp=~/.+(\.txt)/)
  #  	{
  #  		@tempAray=split(/\./,$temp);
  #	  	if(!$locusspecificdata{$tempAray[0]})
  # 	  	{
	#  	  	$locusspecificdata{$tempAray[0]}=$locusspecific."/".$temp;
  #    	}
  #  	}
  #  }
  #  closedir(DIR);
    
  #  foreach my $key (keys %locusspecificdata)
  #  {
  #  	my $OUT=$outpath."/".$key;
  #    &calmethyentrop_locusspecific($locusspecificdata{$key},$OUT);
  #  }
	#}
	if($filetype=~/sf/i)
	{
		&calmethyentrop_locusspecific($locusspecific,$outpath);
	}
	else
	{
		print "Error Number 05.\n";
    print "Unexpected Bug. \n";
    exit(0);
	}
}
else
{
	print "Error Number 05.\n";
  print "Unexpected Bug. \n";
  exit(0);
}

sub calmethynentropy
{
  my ($bismarkdata, $genomepositiondata,$outpath)=@_;
  my $len = 3; ##
  my $reads = 8;##
  my $genomeposition = $genomepositiondata;
  my $temp = "";
  my %genomicCGpositions =();
  my $chrs="";
  my @tempArray = ();
  
  if(!-d $outpath)
	{
		mkdir $outpath;
	}
  
  open(GENOMEPOSITION,  "<$genomeposition" ) or die "\nError Number 02\nCan not open file $genomeposition\n";
  while ($temp=<GENOMEPOSITION>) 
  {
    chomp $temp;
     $temp =~ s/[\r\n]//g;
     $genomicCGpositions{$temp} = $temp;
     $genomicCGpositions{$temp+1} = $temp;
  }
  close(GENOMEPOSITION);
  
  @tempArray = split(/\\/,$bismarkdata);
  @tempArray = split(/\//,$tempArray[scalar(@tempArray)-1]);
  @tempArray = split(/\./,$tempArray[scalar(@tempArray)-1]);
  $chrs = $tempArray[0];
  
  my $inFile = $bismarkdata;
  open(INFILE,  "<$inFile" ) or die "\nError Number 02\nCan not open file $inFile\n";
  my $outFile = $outpath."/".$chrs."_"."readsWith4orMoreCG_CGPOsition_Methycall.txt";
  open (OUTFILE, ">$outFile") or die "\nError Number 02\ncan not open $outFile:$!\n";
  #my $outFile1 = $outpath."/".$chrs."_"."methylationCallWith4orMore_Statistics.xls";
  #open (OUTFILE1, ">$outFile1") or die "\nError Number 02\ncan not open $outFile1:$!\n";

  $temp = "";
  @tempArray = ();
  my $CGcount = 0;
  my %CGcount_his = (); 
  my $seqID = "";
  my @cgIDList = ();
  my $str ="";

  $temp=<INFILE>;
  while ($temp=<INFILE>) 
  {
    chomp $temp;
    $temp =~ s/[\r\n]//g;
    @tempArray = split(/\t/,$temp); 
    my $tempSeqID = $tempArray[0];
    $str = $tempArray[4]; 
    $str =~ s/Z/1/;
    $str =~ s/z/0/;
    if($seqID eq $tempSeqID)
    { 
  	  $CGcount++; 
	    my $tempCGIDmethycalls = join ":",$genomicCGpositions{$tempArray[3]},$str; 
      push(@cgIDList, $tempCGIDmethycalls);
    }
    elsif($seqID ne "")
    {
      if($CGcount >= $len)
      {
	      my @cgIDListsorted  = sort{$a <=>$b} @cgIDList;
	      print OUTFILE "$seqID\t@cgIDListsorted\n";
      }
	    $CGcount_his{$CGcount}++; 
	    @cgIDList = ();
      my $tempCGIDmethycalls = join ":",$genomicCGpositions{$tempArray[3]},$str;
      push(@cgIDList, $tempCGIDmethycalls);                             
	    $seqID = $tempSeqID; 
  	  $CGcount =1;
    }
    else
    {
	    @cgIDList = ();
      my $tempCGIDmethycalls = join ":",$genomicCGpositions{$tempArray[3]},$str;
      push(@cgIDList, $tempCGIDmethycalls);    
    	$seqID = $tempSeqID;
    	$CGcount = 1;
    }
  }
  if($CGcount >= $len)
  {
	  my @cgIDListsorted  = sort{$a <=>$b} @cgIDList;
	  print OUTFILE "$seqID\t@cgIDListsorted\n";
  }
  $CGcount_his{$CGcount}++;
  close(INFILE);
  close(OUTFILE);

  #print OUTFILE1 "4orMoreCGsInReads\t#Count\n";
  #foreach my $key ( sort {$a <=> $b } keys %CGcount_his)
  #{
  #  print OUTFILE1 "$key\t$CGcount_his{$key}\n";
  #}
  #close(OUTFILE1);
  
  my $inFile = $outpath."/".$chrs."_"."readsWith4orMoreCG_CGPOsition_Methycall.txt";
  open( INFILE,  "<$inFile" ) or die "\nError Number 02\nCan not open file $inFile\n";
  #my $outFile = $outpath."/".$chrs."_"."Reads_4orMoreCG_readsCount_16orMoreReads.txt";
  #open (OUTFILE, ">$outFile") or die "\nError Number 02\ncan not open $outFile:$!\n";
  #my $outFile1 = $outpath."/".$chrs."_"."Reads_4orMoreCG_readsCount.txt";
  #open (OUTFILE1, ">$outFile1") or die "\nError Number 02\ncan not open $outFile1:$!\n";
  #my $outFile2 = $outpath."/".$chrs."_"."ReadsCountInSegments_Statistics.xls";
  #open (OUTFILE2, ">$outFile2") or die "\nError Number 02\ncan not open $outFile2:$!\n";
  my $outFile3 = $outpath."/".$chrs."_3CG_16orMoreReads_segments.txt";
  open (OUTFILE_pattern, ">$outFile3") or die "\nError Number 02\ncan not open $outFile3:$!\n";

  $temp = "";
  @tempArray = ();
  @cgIDList = ();
  my %readCount_his = ();
  my $cg4ID = "";
  my %cgIDList_methcall=();
  
  #$temp=<INFILE>;
  while ($temp=<INFILE>) 
  {
    chomp $temp;
    $temp =~ s/[\r\n]//g;
    @tempArray = split(/\t/,$temp);
    @cgIDList = split(/\s+/,$tempArray[1]);
    my $cgIDmethcall = "";
    my $numberOF4CGString = scalar(@cgIDList) -$len;
    for(my $i=0; $i<=$numberOF4CGString; $i++)
    {
      my $a = $i+0;
      my $b = $i+1;
      my $c = $i+2;
      #my $d = $i+3;
      my @cgID1 = split(/:/,$cgIDList[$a]);
      my $cgID1site = $cgID1[0];
      my $cgID1methcall = $cgID1[1];
      my @cgID2 = split(/:/,$cgIDList[$b]);
      my $cgID2site = $cgID2[0];
      my $cgID2methcall = $cgID2[1];
      my @cgID3 = split(/:/,$cgIDList[$c]);
      my $cgID3site = $cgID3[0];
      my $cgID3methcall = $cgID3[1];
      #my @cgID4 = split(/:/,$cgIDList[$d]);
      #my $cgID4site = $cgID4[0];
      #my $cgID4methcall = $cgID4[1];
    
      $cg4ID = $cgID1site."_";
      $cg4ID .= $cgID2site."_";
      $cg4ID .= $cgID3site."_";
      #$cg4ID .= $cgID4site;#."_";
      $cgIDmethcall .= $cgID1methcall;
      $cgIDmethcall .= $cgID2methcall;
      $cgIDmethcall .= $cgID3methcall;
      #$cgIDmethcall .= $cgID4methcall;
    
      $readCount_his{$cg4ID}++;
      if($cgIDList_methcall{$cg4ID})
      {
	      $cgIDList_methcall{$cg4ID} .= "\t".$cgIDmethcall; 
	    }
	    else
	    {
	      $cgIDList_methcall{$cg4ID} = $cgIDmethcall; 
	    }
      $cgIDmethcall= "";
    }
  }
  close(INFILE);
  
  my @cg4IDs_key = ();
  my $keycount=0;
  my %readHis = ();
  
  #print OUTFILE1 "SegmentsID\t#Reads\n";
  #print OUTFILE "SegmentsID\t#Reads(16orMoreReads)\n";
  foreach my $key ( sort {$a <=> $b } keys %readCount_his)
  {
     #print OUTFILE1 "$key\t$readCount_his{$key}\n";
     $temp=$key;
	   if($readCount_his{$key} >= $reads)
	   {
        #print OUTFILE "$key\t$readCount_his{$key}\n";
        $cg4IDs_key[$keycount] = $key;
        $keycount++;        	
	   }
	   $readHis{$readCount_his{$key}}++;
  }

  print OUTFILE_pattern "Title: the Segments Extracted from Methylation Data\n";
  foreach my $key (sort {$a <=> $b } @cg4IDs_key)
  {
	  print OUTFILE_pattern "\n>$key\n";
    my @patternArray = split(/\t/, $cgIDList_methcall{$key});
    for (my $i=0;$i<=$#patternArray; $i++) 
    {
      print OUTFILE_pattern "$patternArray[$i]\n";
    }
  }
  #close ( OUTFILE );
  #close ( OUTFILE1 );
  close (OUTFILE_pattern);

  #print OUTFILE2 "ReadsNumberInSegments\t#Count\n";
  #foreach my $key ( sort {$a <=> $b } keys %readHis)
  #{
  #  print OUTFILE2 "$key\t$readHis{$key}\n"; 
  #}
  #close ( OUTFILE2 );
  
  $temp = ""; 
  @tempArray = (); 
  my $fragment =0; 
  my $regionID = ""; 
  my $fragmentID = ""; 
  my %readPattern = ();

  my $inFile = $outpath."/".$chrs."_3CG_16orMoreReads_segments.txt";;
  open( INFILE,  "<$inFile" ) or die "\nError Number 02\nCan not open file $inFile\n";

  $temp = <INFILE>;
  while($temp = <INFILE>)
  {	
    chomp($temp);
    if($temp =~ />(.*)/)
    {
	    $regionID = $1;	
    }
    elsif($temp ne "")
    {
	    my $mp = $temp;
	    $fragmentID = $regionID;
	    if($mp !~ /2/)
	    {
	      if($readPattern{$fragmentID})
	      {
	        $readPattern{$fragmentID} .= "\t".$mp; 
	      }
	      else
	      {
	        $readPattern{$fragmentID} = $mp; 
	      }
	    }
    }
  }
  close INFILE;

  my $outFile1 = $outpath."/".$chrs."_"."3CG.txt";
  open(OUTFILE1,  ">$outFile1" ) or die "\nError Number 02\nCan not open file $outFile1\n";
  #my $outFile = $outpath."/".$chrs."_"."Chrs_4CG_withMP.txt";
  #open(OUTFILE,  ">$outFile" ) or die "\nError Number 02\nCan not open file $outFile\n";

  print OUTFILE1 "SegmentsID(CG4Postion)\t#Reads\tmethylationLevel\tmethylationEntropy\tmethylatedCG\ttotalCG\n";
  #print OUTFILE "SegmentsID(CG4Postion)\t#ReadsNumber\tmethylationLevel\tmethylationEntropy\tmethylatedCG\ttotalCG\tmethylationPattern\n";
  foreach my $CG4Postion ( sort {$a <=> $b } keys %readPattern)
  {
    my @patternArray = split(/\t/, $readPattern{$CG4Postion});
    my $cgNumber = length($patternArray[0]);
    my %is_newPattern = ();
    my $readNumber = 0;
    my $sumCGMethylated = 0;
    my $sumCGCount = 0;
    my $methylationLevel = 0;
    my $methylationEntropy = 0;

    foreach my $pattern(@patternArray)
    {
      if($is_newPattern{$pattern})
      {
         $is_newPattern{$pattern} ++;
      }
      else
      {
         $is_newPattern{$pattern} = 1;
      }
      $readNumber ++;
    }

    if($readNumber >= $reads)
    {
      #print OUTFILE "$CG4Postion\t$readNumber";
      print OUTFILE1 "$CG4Postion\t$readNumber";
      foreach my $key ( sort {$a <=> $b } keys %is_newPattern)
      {
        my $patternReadCount = $is_newPattern{$key};
        $cgNumber = length($key);
        $sumCGCount += $cgNumber*$patternReadCount;
        $key =~ s/0//g;
        my $cgMethylated = length($key);
        $sumCGMethylated += $cgMethylated*$patternReadCount;
      }
      $methylationLevel = $sumCGMethylated/$sumCGCount;
      foreach my $key ( sort {$a <=> $b } keys %is_newPattern)
      {
        $methylationEntropy += -($is_newPattern{$key}/$readNumber)*log10($is_newPattern{$key}/$readNumber);
      }
      $methylationEntropy *= 3.321928095;
      $methylationEntropy /= $len;
      #printf OUTFILE "\t%.3f\t%.3f\t",$methylationLevel,$methylationEntropy;
      #print OUTFILE "$sumCGMethylated\t$sumCGCount";
      printf OUTFILE1 "\t%.3f\t%.3f\t",$methylationLevel,$methylationEntropy;
      print OUTFILE1 "$sumCGMethylated\t$sumCGCount\n";
      #foreach my $key ( sort {$a <=> $b } keys %is_newPattern)
      #{
      #  print OUTFILE "\t$key\t$is_newPattern{$key}";
      #}
      #print OUTFILE "\n";
    }
  }
  #close(OUTFILE);
  close(OUTFILE1);
  unlink $outpath."/".$chrs."_"."readsWith4orMoreCG_CGPOsition_Methycall.txt";
}

sub calmethyentrop_locusspecific
{
	my ($locusspecific, $outpath)=@_;
  my $temp = ""; 
  my @tempArray = (); 
  my $fragment =0; 
  my $regionID = ""; 
  my $fragmentID = ""; 
  my %readPattern = ();
  my $inFile = $locusspecific;
  my $samplename="";
  
  if(!-d $outpath)
	{
		mkdir $outpath;
	}
	
  @tempArray = split(/\\/,$locusspecific);
  @tempArray = split(/\//,$tempArray[scalar(@tempArray)-1]);
  @tempArray = split(/\./,$tempArray[scalar(@tempArray)-1]);
  $samplename = $tempArray[0];

  open( INFILE,  "<$inFile" ) or die "Error Number 02\nCan not open file $inFile\n";
  while($temp = <INFILE>)
  {	
    chomp($temp);
    if($temp =~ />(.*)/)
    {
	    $regionID = $1;	
    }
    elsif($temp ne "")
    {
      my $cgCount = length($temp);
      $fragment = 0;
      for(my $i = 0; $i <($cgCount-3); $i++)
      {
	      my $mp = substr($temp, $i, 4);
      	$fragmentID = $fragment."_".$regionID;
	      $fragment ++;
	      if($mp !~ /2/)
	      {
	        if($readPattern{$fragmentID})
	        {
	          $readPattern{$fragmentID} .= "\t".$mp; 
	        }
	        else
	        {
	          $readPattern{$fragmentID} = $mp; 
	        }
	      }
      }
    }
  }
  close INFILE;

  my $outFile = $outpath."/".$samplename."_4CG_16orMoreReads"."_segments.txt";
  open(OUTFILE,  ">$outFile" ) or die "Error Number 02\nCan not open file $outFile\n";
  my $outFile1 = $outpath."/".$samplename."_4CG.txt";
  open(OUTFILE1,  ">$outFile1" ) or die "Error Number 02\nCan not open file $outFile1\n";
  print OUTFILE "Title: the Segments Extracted from Methylation Data.\n";
  print OUTFILE1 "LocusID\t#Number\tLocus\t#Reads\tMethylationLevel\tMethylationEntropy\tmethylatedCG\ttotalCG\n";
  foreach my $CG4Postion ( sort {$a cmp $b } keys %readPattern)
  {
    my @patternArray = split(/\t/, $readPattern{$CG4Postion});
    my $cgNumber = length($patternArray[0]);
    my %is_newPattern = ();
    my $readNumber = 0;
    my $sumCGMethylated = 0;
    my $sumCGCount = 0;
    my $methylationLevel = 0;
    my $methylationEntropy = 0;
    

    foreach my $pattern(@patternArray)
    {
      if($is_newPattern{$pattern})
      {
        $is_newPattern{$pattern} ++;
      }
      else
      {
        $is_newPattern{$pattern} = 1;
      }
      $readNumber ++;
    }

    if($readNumber >= 16)
    {
      if($CG4Postion =~ /(.*)_(.*)/)
      {
        print OUTFILE1 "$CG4Postion\t$1\t$2\t$readNumber";
      }
      foreach my $key ( sort {$a <=> $b } keys %is_newPattern)
      {
        my $patternReadCount = $is_newPattern{$key};
        $cgNumber = length($key);
        $sumCGCount += $cgNumber*$patternReadCount;
        $key =~ s/0//;
        my $cgMethylated = length($key);
        $sumCGMethylated += $cgMethylated*$patternReadCount;
      }
      $methylationLevel = $sumCGMethylated/$sumCGCount;
      foreach my $key ( sort {$a <=> $b } keys %is_newPattern)
      {
        $methylationEntropy += -($is_newPattern{$key}/$readNumber)*log10($is_newPattern{$key}/$readNumber);
      }
      $methylationEntropy *= 3.321928095;
      $methylationEntropy /= 4;
      
      @tempArray = split(/\t/,$readPattern{$CG4Postion});
      $readPattern{$CG4Postion} = join "\n", @tempArray;
      print OUTFILE "\n>$CG4Postion\n";
      print OUTFILE "$readPattern{$CG4Postion}\n";
      printf OUTFILE1 "\t%.3f\t%.3f",$methylationLevel,$methylationEntropy;
      print OUTFILE1 "\t$sumCGMethylated\t$sumCGCount\n";
    }
  }
  close(OUTFILE);
  close(OUTFILE1);
}

sub log10 
{
  my $n = shift;
  return log($n)/log(10);
}
