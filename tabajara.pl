#!/usr/bin/perl

# tabajara.pl - A tool for rational design of profile HMMs
# v1.02: Adding a routine to calculate a minimum cutoff score for profile HMMs constructed from full-length MSA.
#        Change in model validation routine: use of domain score value to verify if the constructed model is valid.
# v1.01: Removing unecessary comments.

use strict;
use Getopt::Long;
use File::Basename;
use constant  PSEUDOCOUNT => 0.0000001;

# turn on buffer autoflush
$| = 1;

# Capture the commands given to Tabajara in the command line 
my $script_command = $0;
foreach (@ARGV) {
    $script_command .= /\s/ ?   " \"" . $_ . "\""
                    :           " "   . $_;
}

#variables
my $version = "1.02";
my $last_update = "2021-08-09";
my $version_option;
my $window_size_score;
my $win_lam;
my $outfile_name;
my @bg_distribution;
my $distribution;
my $scoring_function;
my $use_seq_weights;
my $gap_cutoff;
my $use_gap_penalty;
my $seq_specific_output;
my $normalize = "yes";
my $input_file;
my $help_print;
my $help;
my $output = "output_dir";
my $len;
my @seq_weights;
my @amino_acids = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-');
my @iupac_alphabet = ("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-");
my @nucleotides = ('A', 'C', 'G', 'T');
my $window_size;
my $porc_threshold;
my $threshold;
my $min_size_block;
my @values;
my @valid_blocks;
my @valid_windows;
my @final_blocks;
my $i;
my $j;
my $start;
my $end;
my $count;
my $start_w;
my $end_w;
my $aux;
my @aux2;
my $len_values;
my $interval;
my $help_string;
my $help;
my $log;
my $gap_perc;
my $max_block_size;
my $amount = 2;
my $measure;
my $category;
my $remove_redundance = "no"; 
my $hmmbuild;
my $cat;
my $window_gaps = "no";
my $threshold_validation;
my $amount_detection = 80;
my $category_percentage = 80;
my $pontuation = 1;
my $clean = "yes";
my $cutoff_score = "no";
my $more_blocks = "no";
my $conf;
my %config;
my $full_length = "no";

GetOptions ("i|input_file=s"			=>\$input_file,
	    "a=s"				=>\$seq_specific_output,
	    "o|output=s"			=>\$output,
	    "h|help"				=>\$help,
	    "l=s"				=>\$use_seq_weights,
	    "gs|gap_sequence=s"			=>\$gap_perc,
	    "gp=s"				=>\$use_gap_penalty,
	    "sa|sequence_amount=s"		=>\$amount,
	    "fl|full_length=s"                    =>\$full_length,
	    "gc|gap_cutoff=s"			=>\$gap_cutoff,
	    "n=s"				=>\$normalize,
	    "s=s"				=>\$scoring_function,
	    "w|window_size=i"       		=>\$window_size,
            "p|percentage=f"       		=>\$porc_threshold,
            "t|threshold=f"       		=>\$threshold,
            "b|minimum_block_size=i"   		=>\$min_size_block,
            "md|maximum_distance=s" 		=>\$interval,
	    "mb|maximum_block_size=s"		=>\$max_block_size,
	    "m|method=s"			=>\$measure,
	    "c|category=s"			=>\$cat,
	    "wg|discard_windows=s"     		=>\$window_gaps,
	    "v|version"				=>\$version_option,
	    "di|discard_sequences=s"		=>\$remove_redundance,
	    "pt|maximum_percentage_ratio=f"	=>\$threshold_validation,
	    "pd|minimum_detection_rate=f"	=>\$amount_detection,
	    "pc|minimum_category_rate=f"	=>\$category_percentage,
	    "sv|score_value=f"			=>\$pontuation,
	    "hb|hmmbuild_parameters=s"		=>\$hmmbuild,
	    "sr|search_more_regions=s"		=> \$more_blocks,
	    "clean=s"				=>\$clean,
	    "conf=s"    			=>\$conf,
	    "cs|cutoff_score=s"			=>\$cutoff_score);

$help_print = " 
TABAJARA - a tool for rational design of profile HMMs - $version ($last_update)
(c) 2020. Liliane S. Oliveira & Arthur Gruber
Usage: tabajara.pl -i <input_file> -o <output_directory>

Mandatory parameters:
-i|input_file <file name>	 	Input file (multiple alignment in FASTA or CLustal formats).
-m|method <c|d>                         Method to select blocks. Options:
                                         c - Conservation
                                         d - Discrimination

NOTE: When using -m d, the parameter -c becomes mandatory, otherwise this parameter is optional:

-c|category <string>                    Name of major category to be analysed.

IMPORTANT:

When using -fl no (default), the following parameters are also mandatory:

-b|minimum_block_size <integer>         Minimum block size (must be ≥ w).
-t|threshold <decimal>    		Threshold of the alignment position for block extraction  (valid values: 0 to 1).
-p|percentage <integer>    		Percentage of positions in sliding window with score ≥ t.
-w|window_size <integer>   		Window size for block extraction.


Optional parameters:
-conf           			Configuration file
-clean <yes|no>				Remove all intermediate files during the execution of Tabajara (default = yes).
-cs|cutoff_score <yes|no>               Insert cutoff scores in the profile HMMs (default = no).
-di|discard_sequences <yes|no>    	Discard identical sequences (default = no).
-fl|full_length <yes|no>                Use full-length sequence for model construction (default = no).
-gc|gap_cutoff <integer>   		Gap cutoff. Do not score columns that contain more than gap cutoff fraction gaps (default = 30).
-gs|gap_sequence <decimal>   		Percentage of allowed gaps in each sequence.
-hb|hmmbuild_parameters <string>    	Any set of hmmbuild's valid parameters can be entered under double quotes. 
					    Example: -hb \"--wblosum --wid 0.8\". Default: no parameters.
-h|help              			Help.
-mb|maximum_block_size <integer>	Maximum block size.
-md|maximum_distance <integer>   	Maximum accepted distance for two alignment blocks to be joined as a single block.
-o|output              			Output directory.
-pc|minimum_category_rate		Minimum percentage of categories in the training set that must meet the criteria defined by the parameter pd for profile HMMs built in conservation mode (default = N).
-pd|minimum_detection_rate <decimal>   	Minimum detection rate (in percentage) of training set (MSA) sequences by the constructed HMM. 
					This is equivalent to the minimum accepted sensitivity of the model (default = 80).
-pt|maximum_percentage_ratio <decimal>  This parameter applies for the validation of HMMs designed for group discrimination. Given two groups, 
					A (selected category) and B (remaining sequences), an hmmsearch is performed with the constructed model 
					against the training set (MSA) sequences. The highest score obtained by group B sequences must be lower 
					than n % (defined by parameter -pt) of the lowest score obtained by group A sequences, otherwise the model 
					is not accepted (default = 80).
-sa|sequence_amount <integer>           Minimum sequence amount of MSA to construct hmm (default = 2).
-sr|search_more_regions <yes|no>	This parameter is used to search for more subregions in a previously selected region (default = no). 
-sv|score_value <decimal>		This parameter is used to calculate the minimum score for HMMs designed by Tabajara. The minimum score will 
					be used as suggested score if calculated score is lower than this score. This minimum score is calculated by 
					mutiplying the value of -sv by HMM length (default = 1).
-v|version				Version.
-wg|discard_windows <yes|no>    	Discard sliding windows presenting gaps (default = no).
\n";

if($help){
    die $help_print;
}

if($version_option){
    die "Version $version\nLast update: $last_update\n";
}

#
# if configuration file is not defined, check if the mandatory arguments are defined 

if (not defined $conf) {
    if(!$input_file){
    	die "ERROR: Missing mandatory argument -i.\n$help_print\n";
    }

    if(!$measure){
    	die "ERROR: Missing mandatory argument -m.\n$help_print\n";
    }

    if(lc($full_length) eq "no"){
    	if(!$porc_threshold){
    	    die "ERROR: Missing mandatory argument -p.\n$help_print\n";
    	}

    	if($threshold == 0){
	    $threshold = 0;
    	}
    	elsif(!$threshold){
    	    die "ERROR: Missing mandatory argument -t.\n$help_print\n";
    	}

    	if(!$min_size_block){
            die "ERROR: Missing mandatory argument -b.\n$help_print\n";
    	}

    	if(!$window_size){
            die "ERROR: Missing mandatory argument -w.\n$help_print\n";
    	}
    }
}

# If configuration file is defined, check if all mandatory arguments are defined 

my @strs = ();
if ($conf) {
    print STDERR "Configuration file specified, command line options will be overridden.\n";

    open(CONFIG, "< $conf") or die("ERROR: Problem opening configuration file $conf: $!\n");

    my $configLine;
    while ($configLine = <CONFIG>) {
        $configLine =~ s/^\s+//;
        $configLine =~ s/\s+\Z//;
        $configLine =~ s/\s+\=/\=/;
        $configLine =~ s/\=\s+/\=/;
        if ($configLine =~ /^\#/ || !($configLine =~ /(.)\=(.)/)) {
            next;
        }
        chomp $configLine;

        if ($configLine =~ m/(.+?)=(.+)/) {
            $config{$1} = $2;
        }
    }
    close(CONFIG);
   
    if (defined($config{"fl"}) and defined($config{"full_length"})){
        print "Parameters -fl and -full_length defined. Tabajara will consider the value of -fl parameter\n";
        $full_length = $config{"fl"};
    }
    else{
        if (defined($config{"fl"})){
            $full_length = $config{"fl"};
        }
        elsif(defined($config{"full_length"})){
            $full_length = $config{"full_length"};
        }
    }
    
    #
    # Mandatory arguments
    #
    my $missingArgument = 0;
    if(lc($full_length) eq "no"){
    	if (!$config{"b"} and !$config{"minimum_block_size"}) {
            $missingArgument = 1;
            print "Missing mandatory configuration argument: Minimum block size (-b).\n";
    	}
    	else {
	    if($config{"b"}){
            	$min_size_block = $config{"b"};
	    }
	    elsif($config{"minimum_block_size"}){
            	$min_size_block = $config{"minimum_block_size"};
	    }
	}
    }

    if (!$config{"i"} and !$config{"input_file"}) {
        $missingArgument = 1;
        print "Missing mandatory configuration argument: input file (-i).\n";
    }
    else {
	if($config{"i"}){
            $input_file = $config{"i"};
	}
	elsif($config{"input_file"}){
            $input_file = $config{"input_file"};
        }
    }

    if (!$config{"m"} and !$config{"method"}) {
        $missingArgument = 1;
        print "Missing mandatory configuration argument: method (-m).\n";
    }
    else {
	if($config{"m"}){
            $measure = $config{"m"};
	}
	elsif($config{"method"}){
	    $measure = $config{"method"};
	}
    }
    if(lc($full_length) eq "no"){
    	if (!$config{"t"} and !$config{"threshold"}) {
            $missingArgument = 1;
            print "Missing mandatory configuration argument: threshold (-t).\n";
    	}
    	else {
	    if($config{"t"}){
            	$threshold = $config{"t"};
	    }
	    elsif($config{"threshold"}){
            	$threshold = $config{"threshold"};
            }
    	}

    	if (!$config{"percentage"} and !$config{"p"}) {
            $missingArgument = 1;
            print "Missing mandatory configuration argument: percentage (-p).\n";
    	}
    	else {
	    if($config{"p"}){
            	$porc_threshold = $config{"p"};
	    }
	    elsif($config{"percentage"}){
	    	$porc_threshold = $config{"percentage"};
	    }
    	}

    	if (!$config{"w"} and !$config{"window_size"}) {
            $missingArgument = 1;
            print "Missing mandatory configuration argument: window size (-w).\n";
    	}
    	else {
	    if($config{"w"}){
           	$window_size = $config{"w"};
	    }
	    elsif($config{"window_size"}){
	    	$window_size = $config{"window_size"};
	    }
    	}
    }
    if ($missingArgument) {
        die "\nERROR: Cannot run tabajara.pl, mandatory configuration argument(s) missing (see bellow).\n\n$help_print\n";
    }

    #
    # Optional parameters
    #

    my $max_ent = -5*((1/5)*log2(1/5));
    if (defined($config{"c"}) and defined($config{"category"})){
	print "Parameters -c and -category defined. Tabajara will consider the value of -c parameter\n";
        $cat = $config{"c"};
    }
    else{
	if (defined($config{"c"})){
	    $cat = $config{"c"};
	}
	elsif(defined($config{"category"})){
            $cat = $config{"category"};
        }
    }

    if (defined($config{"clean"})){
        $clean = $config{"clean"};
    }
    
    if (defined($config{"di"}) and defined($config{"discard_sequences"})){
	print "Parameters -di and -discard_sequences defined. Tabajara will consider the value of -di parameter\n";
        $remove_redundance = $config{"di"};
    }
    else{
	if(defined($config{"di"})){
	    $remove_redundance = $config{"di"};
	}
	elsif(defined($config{"discard_sequences"})){
            $remove_redundance = $config{"discard_sequences"};
        }
    }

    if(lc($full_length) eq "no"){
    	if (defined($config{"gc"}) and defined($config{"gap_cutoff"})){
	    print "Parameters -gc and -gap_cutoff defined. Tabajara will consider the value of -gc parameter\n";
            $gap_cutoff = $config{"gc"};
    	}
    	else{
	    if (defined($config{"gc"})){
	    	$gap_cutoff = $config{"gc"};
	    }
	    elsif(defined($config{"gap_cutoff"})){
	    	$gap_cutoff = $config{"gap_cutoff"};
	    }
    	}

    	if (defined($config{"gs"}) and defined($config{"gap_sequence"})){
	    print "Parameters -gs and -gap_sequence defined. Tabajara will consider the value of -gs parameter\n";
            $gap_perc = $config{"gs"};
    	}
    	else{
	    if (defined($config{"gs"})){
            	$gap_perc = $config{"gs"};
            }
	    elsif(defined($config{"gap_sequence"})){
            	$gap_perc = $config{"gap_sequence"};
            }
    	}

	if (defined($config{"mb"}) and defined($config{"maximum_block_size"})){
            print "Parameters -mb and -maximum_block_size defined. Tabajara will consider the value of -mb parameter\n";
            $max_block_size = $config{"mb"};
    	}
    	else{
            if(defined($config{"mb"})){
            	$max_block_size = $config{"mb"};
            }
            elsif(defined($config{"maximum_block_size"})){
            	$max_block_size = $config{"maximum_block_size"};
            }
    	}

    	if (defined($config{"md"}) and defined($config{"maximum_distance"})){
           print "Parameters -md and -maximum_distance defined. Tabajara will consider the value of -md parameter\n";
           $interval = $config{"md"};
    	}
    	else{
            if(defined($config{"md"})){
            	$interval = $config{"md"};
            }
            elsif(defined($config{"maximum_distance"})){
            	$interval = $config{"maximum_distance"};
            }
    	}
	if (defined($config{"sa"}) and defined($config{"sequence_amount"})){
            print "Parameters -sa and -sequence_amount defined. Tabajara will consider the value of -sa parameter\n";
            $amount = $config{"sa"};
    	}
    	else{
            if (defined($config{"sa"})){
            	$amount = $config{"sa"};
            }
            elsif(defined($config{"sequence_amount"})){
            	$amount = $config{"sequence_amount"};
            }
    	}
	if (defined($config{"sv"}) and defined($config{"score_value"})){
            print "Parameters -sc and -score_value defined. Tabajara will consider the value of -sv parameter\n";
            $pontuation = $config{"sv"};
    	}
    	else{
            if (defined($config{"sv"})){
            	$pontuation = $config{"sv"};
            }
            elsif(defined($config{"score_value"})){
            	$pontuation = $config{"score_value"};
            }
    	}

    	if (defined($config{"wg"}) and defined($config{"discard_windows"})){
            print "Parameters -wg and -discard_windows defined. Tabajara will consider the value of -wg parameter\n";
            $window_gaps = $config{"wg"};
    	}
    	else{
            if (defined($config{"wg"})){
            	$window_gaps = $config{"wg"};
            }
            elsif(defined($config{"discard_windows"})){
            	$window_gaps = $config{"discard_windows"};
            }
    	}

        if (defined($config{"sr"}) and defined($config{"search_more_regions"})){
            print "Parameters -sr and -search_more_regions defined. Tabajara will consider the value of -sr parameter\n";
            $more_blocks = $config{"sr"};
    	}
    	else{
            if (defined($config{"sr"})){
            	$more_blocks = $config{"sr"};
            }
            elsif(defined($config{"search_more_regions"})){
            	$more_blocks = $config{"seach_more_regions"};
            }
    	}
    }

    if (defined($config{"hb"}) and defined($config{"hmmbuild_parameters"})){
	print "Parameters -hb and -hmmbuild_parameters defined. Tabajara will consider the value of -hb parameter\n";
        $hmmbuild = $config{"hb"};
    }
    else{
	if (defined($config{"hb"})){
	    $hmmbuild = $config{"hb"};
	}
	elsif(defined($config{"hmmbuild_parameters"})){
	    $hmmbuild = $config{"hmmbuild_parameters"};
	}
    }

    if (defined($config{"o"}) and defined($config{"output"})){
	print "Parameters -o and -output defined. Tabajara will consider the value of -o parameter\n";
        $output = $config{"o"};
    }
    else{
	if(defined($config{"o"})){
	    $output = $config{"o"};
	}
	elsif(defined($config{"output"})){
	    $output = $config{"output"};
	}
    }

    if (defined($config{"pd"}) and defined($config{"minimum_detection_rate"})){
	print "Parameters -pd and -minimum_detection_rate defined. Tabajara will consider the value of -pd parameter\n";
        $amount_detection = $config{"pd"};
    }
    else{
	if(defined($config{"pd"})){
            $amount_detection = $config{"pd"};
        }
        elsif(defined($config{"minimum_detection_rate"})){
            $amount_detection = $config{"minimum_detection_rate"};
        }
    }

   if (defined($config{"pc"}) and defined($config{"minimum_category_rate"})){
        print "Parameters -pc and -minimum_category_rate defined. Tabajara will consider the value of -pc parameter\n";
        $category_percentage = $config{"pc"};
    }
    else{
        if(defined($config{"pc"})){
            $category_percentage = $config{"pc"};
        }
        elsif(defined($config{"minimum_category_rate"})){
            $category_percentage = $config{"minimum_category_rate"};
        }
    }

    if (defined($config{"pt"}) and defined($config{"maximum_percentage_ratio"})){
	print "Parameters -pt and -maximum_percentage_ratio defined. Tabajara will consider the value of -pt parameter\n";
        $threshold_validation = $config{"pt"};
    }
    else{
	if (defined($config{"pt"})){
	    $threshold_validation = $config{"pt"};
	}
	elsif(defined($config{"maximum_percentage_ratio"})){
	    $threshold_validation = $config{"maximum_percentage_ratio"}
	}
    }

    if (defined($config{"cs"}) and defined($config{"cutoff_score"})){
	print "Parameters -sc and -score_cutoff defined. Tabajara will consider the value of -sc parameter\n";
        $cutoff_score = $config{"cs"};
    }
    else{
	if(defined($config{"cs"})){
	    $cutoff_score = $config{"cs"};
	}
	elsif(defined($config{"cutoff_score"})){
	    $cutoff_score = $config{"cutoff_score"};
	}
    }
}

if((lc($measure) ne "mi") and (lc($measure) ne "sh") and (lc($measure) ne "d") and (lc($measure) ne "c")){
    die "ERROR: Invalid method! Select mi, sh, d or c!\n$help_print\n";
}

if(lc($full_length) eq "no"){
    if($porc_threshold < 0 or $porc_threshold > 100){
    	 die "ERROR: Invalid value! Valid values: 0 to 100.\n$help_print\n";
    }

    if($threshold > 1 or $threshold < 0){
    	die "ERROR: Invalid value! Valid values: 0 to 1 .\n$help_print\n";
    }

    if(!$min_size_block){
    	die "ERROR: Missing mandatory argument -b.\n$help_print\n";
    }

    if($window_size <= 0){
    	die "Invalid value! Window size must be bigger than zero!\n$help_print\n"
    }

    if($window_size > $min_size_block){
    	print "Minimum block size (-b) must be higher or equal to window size (-w)! -b = $window_size\n ";
    	$min_size_block = $window_size; 
    }

    if(lc($more_blocks) ne "yes" and lc($more_blocks) ne "no"){
    	print "Search more regions (-sr) must be yes or no| Ignoring informed value (-sr = no)\n";
    	$more_blocks = "no";
    }
}
else{
    if($porc_threshold){
	push @strs, "-p";
    }

    if($threshold){
	push @strs, "-t";
    }

    if($min_size_block){
	push @strs, "-b";
    }

    if($window_size){
	push @strs, "-w";
    }
}

my $file_type = verifyFileType($input_file);

$clean = lc($clean);

if(!($clean eq "yes") and !($clean eq "no")){
    print "The value of -clean parameter must be yes or no! Ignoring informed value (-clean = yes)\n";
    $clean = "yes";
}

$cutoff_score = lc($cutoff_score);
if(!($cutoff_score eq "yes") and !($cutoff_score eq "no")){
    die "Invalid value! The value must be yes or no!\n$help_print\n";
}

$full_length = lc($full_length);
if(!($full_length eq "yes") and  !($full_length eq "no")){
    print "The value of -fl parameter must be yes or no! Ignoring informed value (-fl = no)\n";
    $full_length = "no";
}

if($use_gap_penalty){
    if(lc($use_gap_penalty) eq 'false'){
	$use_gap_penalty = 0;
    }
}
else{
   $use_gap_penalty = 1;
}

if(!$window_size_score){
   $window_size_score = 3;
}

if(!$win_lam){
   $win_lam = 0.5;
}

if(lc($full_length) eq "no"){
   if($interval){
   	if($interval =~ /^[+-]?\d+$/){
   	}
   	else{
	    die "Invalid value! Maximum accepted distance (-md) must be a number \n$help_print\n";
   	}
	
   	if($interval < 0){
	    die "Invalid value! Maximum accepted distance (-md) must be bigger than zero!\n$help_print\n"
   	}
    }

    if($max_block_size){
   	if($max_block_size =~ /^[+-]?\d+$/){
   	}
   	else{
            die "Invalid value! Maximum block size (-mb)  must be a number!\n$help_print\n";
   	}

       	if($max_block_size < 0){
            die "Invalid value! Maximum block size (-mb) must be bigger than zero!\n$help_print\n"
   	}
    }

   if($amount){
   	if($amount =~ /^[+-]?\d+$/){
   	}
   	else{
            die "Invalid value! Minimum sequence amount (-sa)  must be a number!\n$help_print\n";
   	}

       	if($amount < 0){
            die "Invalid value! Minimum sequence amount (-sa) must be bigger than zero!\n$help_print\n"
   	}

    }

    if($gap_cutoff){
   	if($gap_cutoff =~ /^[+-]?\d+$/){
       	}
   	else{
            die "Invalid value! Gap cutoff (-gc) must be a integer number!\n$help_print\n";
   	}

   	if($gap_cutoff < 0){
            die "Invalid value! Gap cutoff (-gc) must be bigger than zero!\n$help_print\n"
   	}
    }

    if($gap_perc){
   	if($gap_perc =~ /^[+-]?\d+$/){
   	}	
   	else{
            die "Invalid value! Gap cutoff in each sequence (-gs) must be a number!\n$help_print\n";
   	}

   	if($gap_perc < 0){
            die "Invalid value! Gap cutoff in each sequence (-gs) must be bigger than zero!\n$help_print\n"
   	}
    }
}

if($remove_redundance){
    if((lc($remove_redundance) ne "yes") and (lc($remove_redundance) ne "no")){
	die "Invalid value! Discard identical sequences (-di)  must be yes or no!\n$help_print\n";
    }
}

if((!defined $threshold_validation) and (lc($measure) ne "c")){
     $threshold_validation = 80;
}
elsif((defined $threshold_validation) and (lc($measure) eq "c")){
    print STDERR "Invalid parameter -pt for method \"conservation\" (-m c). Ignoring parameter -pt\n";
}

if($amount_detection){
   if($amount_detection =~ /^[+-]?\d+$/){
   }
   else{
        die "Invalid value! The percentage (-pd) must be a number!\n$help_print\n";
   }

   if($amount_detection < 0){
        die "Invalid value! -pd must be bigger than zero!\n$help_print\n"
   }
}

if($threshold_validation){
   if($threshold_validation =~ /^[+-]?\d+$/){
   }
   else{
        die "Invalid value! Maximum percentage ratio (-pt) must be a number!\n$help_print\n";
   }

   if($amount < 0){
        die "Invalid value! Maximum percentage ratio (-pt) must be bigger than zero!\n$help_print\n"
   }

}

if($distribution){
    open(FILE, "$distribution");
    while(<FILE>){
	if($_ =~ /#/){}
	else{
	    @bg_distribution = split(" ", $_);
	}
    }
    close(FILE);
}
else{
    if($file_type == 1){#DNA sequences
   	@bg_distribution = (0.25, 0.25, 0.25, 0.25);
    }
    else{#Protein sequences
    	@bg_distribution = (0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072);
    }
}

if(!$scoring_function){
    $scoring_function = 'js_divergence';
}

if(!$use_seq_weights){
    $use_seq_weights = 'true';
}

if(!$gap_cutoff){
    $gap_cutoff = 30;
}

# Creates a matrix with amino acid counts to each alignment column 
my %aa_to_index;
my $count = 0;
foreach my $aa (@amino_acids){
    $aa_to_index{$aa} = $count;
    $count++;
}

# Creates the output directory
my $caracter = substr $output, -1;
if($caracter eq '/'){
    my $l = length($output);
    $output = substr $output, 0, ($l-1);
}
$output = output_dir_name($output);
system "mkdir $output";

$len = length($input_file);
my $aux_name = substr $input_file, (rindex($input_file, "/") + 1), $len;
$outfile_name = substr $aux_name, 0, index($aux_name, ".");

###############################################################################
#                             Begin execution                                 #
###############################################################################

$log = "logfile.txt";
my $log_file = $output."/".$log;
open(LOG, ">$log_file");

my @aux_print = split(" ",$script_command);
$script_command =~ s/$aux_print[0]/tabajara.pl/;
print LOG "Command line:\n$script_command\n\n";

# Check the composition (nucleotide or amino acid) of the fasta file 
print LOG "Sequence type: ";
if($file_type == 1){
    print LOG "DNA\n";
}
else{
    print LOG "Protein\n";
}

# Verifies the method to be used in the analysis (c - Conservation; d - Discrimination)
print LOG "\n-m Method: ";
if(lc($measure) eq "d" ){
    if(!$cat){
	system "rm -rf $output";
        die "ERROR: Missing mandatory argument -d.\n$help_string\n";
    }
    print LOG "Discrimination.\n";
}
elsif(lc($measure) eq "c"){
    print LOG "Conservation.\n";
}
if($cat){
    print LOG "-c Category(ies): $cat.\n";
}

my @categories = ();

# Checks if categories should be used in the analysis.
if($cat){	
    my @aux_cat = split(",", $cat);
    my %hash_cat = ();
    for(my $i = 0; $i < scalar(@aux_cat); ++$i){
	if(!defined $hash_cat{$aux_cat[$i]}){
	    $hash_cat{$aux_cat[$i]} = 1;
	}
    }
    foreach my $key (sort keys %hash_cat){
	push @categories, $key;
    } 
    %hash_cat = ();
}

my @amount_cat = ();

# Adding the parameters to be used in the analysis in the log file.

print LOG "\nParameters:\n\n";

print LOG "Input/output:\n";
print LOG "-i Input file: $input_file\n";
print LOG "-o Output directory: $output\n\n";

if(lc($full_length) eq "no"){
    print LOG "Block extraction:\n";
    print LOG "-w Window size = $window_size\n";
    print LOG "-p Minimum % of positions with score ≥ t = $porc_threshold\n";
    print LOG "-t Threshold value for block extraction = $threshold\n";
    print LOG "-b Minimum block size = $min_size_block\n\n";
}
else{
    if(defined $window_size){
        push @strs, "-w";
    }
    if(defined $porc_threshold){
        push @strs, "-p";
    }
    if(defined $threshold){
        push @strs, "-t";
    }
    if(defined $min_size_block){
        push @strs, "-b";
    }
}

if(lc($full_length) eq "no"){
    print LOG "Other block extraction parameters:\n";
}
else{
    print LOG "Profile HMM construction parameters:\n";
}
print LOG "-di Discard identical sequences = $remove_redundance ";
if($remove_redundance eq "no"){
    print LOG "(default)";
}
print LOG "\n"; 
if(lc($full_length) eq "yes"){
    print LOG "-fl Use full-length sequence for model construction = yes\n"
}

if(lc($full_length) eq "no"){
    print LOG "-fl Use full-length sequence for model construction = no\n";
    print LOG "-gc Percentage of allowed gaps in each column = $gap_cutoff";
    if($gap_cutoff == 30){
    	print LOG " (default)";
    }
    print LOG "\n";
    print LOG "-gs Percentage of allowed gaps in each sequence = ";
    if($gap_perc){
        print LOG "$gap_perc\n";
    }
    else{
        print LOG "not defined\n";
    }
}
else{
    if($gap_cutoff > 30){
	push @strs, "-gc";
    }
    if($gap_perc){
        push @strs, "-gs";
    }	
}


print LOG "-hb Hmmbuild parameters = ";
if(defined $hmmbuild){
    print LOG "$hmmbuild\n";
}
else{
    print LOG "not defined\n";
}

if(lc($full_length) eq "no"){
    print LOG "-mb Maximum block size = ";
    if($max_block_size){
    	print LOG "$max_block_size\n";
    }
    else{
    	print LOG "not defined\n";
    }
    print LOG "-md Maximum interval between blocks = ";
    if($interval){
    	print LOG "$interval\n";
    }
    else{
    	print LOG "not defined\n"
    }
    print LOG "-sa Minimum # of sequences to construct HMM = $amount";
    if($amount == 2){
    	print LOG " (default)";
    }
    print LOG "\n";
    print LOG  "-wg Discard sliding window(s) containing gaps = $window_gaps";
    if($window_gaps eq "no"){
    	print LOG " (default)";
    }
    print LOG "\n";
    print LOG "-sr Search more regions = $more_blocks";
    if(lc($more_blocks) eq "no"){
    	print LOG " (default)";
    }
}
else{
    if($max_block_size){
        push @strs, "-mb";
    }
    if($interval){
        push @strs, "-md";
    }
   if($amount > 2){
        push @strs, "-sa";
   }
   if(lc($more_blocks) eq "yes"){
        push @strs, "-sr";
   }
   if($window_gaps eq "yes"){
        push @strs, "-wg";
   }
}
my $frase = "";
if(scalar(@strs) > 0){
    if(scalar(@strs) == 1){
    	$frase = "The parameter $strs[0] is not applicable to full protein analysis. This parameter will be ignored!";
    }
    else{
    	$frase = "The parameters $strs[0]";
    	my $aux_size = scalar(@strs);
    	for(my $s = 1; $s < $aux_size - 1; ++$s){
            $frase .= ", $strs[$s]";
    	}
    	$frase .= " and $strs[$aux_size-1] are not applicable to full protein analysis. These parameters will be ignored!"
    }
}
print LOG $frase;
print LOG "\n\n";

print LOG "Profile HMM validation:\n";
print LOG "-pd Minimum percentage of detected sequences = $amount_detection\n";
print LOG "-pc Minimum pergentage of categories with at least -pd percent of detected sequences = $category_percentage \n";

if($measure eq "c"){
    print LOG "-pt Does not apply to method Conservation\n";
}
else{
    print LOG "-pt Maximum score percentage of non-chosen group sequences = $threshold_validation";
    if($threshold_validation == 80){
	print LOG " (default)";
    }
    print LOG "\n";
}
print LOG "-sv Value for minimum score calculation = $pontuation";
if($pontuation == 1){
    print LOG " (default)";
}
print LOG "\n";
print LOG "-cs Insert cutoff scores in the profile HMMs = $cutoff_score\n\n";

$len = length($input_file);
my $align_suffix = substr $input_file, (rindex($input_file, ".") + 1), $len;
my @names = ();
my @alignment = ();
my @seq_weigths = ();

# Reads the aligment file

# Identifies the format of the alignment file (fasta or clustal)
my $t = verify_file_type($input_file);
my $aux1;
my $aux2;
my $aux3;

if($t == 1){# fasta file
    ($aux1, $aux2, $aux3) = read_fasta_alignment($input_file);
}
else{# clustal file 
    ($aux1, $aux2, $aux3) = read_clustal_alignment($input_file);
}

my $pos = index($input_file, ".");
my $pos_f = rindex($input_file, "/");
my $input_file_prefix; 
if($pos_f > -1){
    $pos -= ($pos_f+1);
    $input_file_prefix = substr $input_file, ($pos_f+1), ($pos);
}
else{
    $input_file_prefix = substr $input_file, 0, ($pos);
}

@names = @{$aux1};
@alignment = @{$aux2};

if(@names eq ()){
    print LOG "Could not identify $input_file file format. Exiting...\n";
    close(LOG);
    die "Could not identify $input_file file format. Exiting...";
}
if (scalar(@alignment) != scalar(@names) or @alignment eq ()){
    print LOG "Unable to parse alignment.\n";
    close(LOG);
    die "Unable to parse alignment.\n";
}

my $no_red = $output."/".$outfile_name;
print LOG "Pre-processing:\n";

if(lc($remove_redundance) eq "yes"){

    #Removes redundant sequences from the aligment file

    my $count;
    $no_red .="_no_redundancy";
    ($count, $aux1, $aux2, $aux3) = remove_redundancy_MSA(\@{$aux1}, \@{$aux2});
    if($count == 0){
    	print LOG "No identical sequences found.\n";
    }
    else{
    	print LOG "$count identical sequences removed.\n";
    }

    @names = @{$aux1};
    @alignment = @{$aux2};
}

# Checks the number of sequences in aligment file 

my $original_size = scalar(@names);
if($original_size < $amount){
    print LOG "Number of sequences in $input_file is lower than $amount\n";
    close(LOG);
    die "Number of sequences in $input_file is lower than $amount\n";
}
else{
    if($original_size == 1){
    	if (lc($remove_redundance) eq "yes"){
	    print LOG "$input_file contains just 1 sequence. The sequence will be repeated 5 times in the file and parameter -di = no\n";
	    $remove_redundance = "no";
            print STDERR "$input_file contains just 1 sequence. The sequence will be repeated 5 times in the file and parameter -di = no\n";
        }
    	for(my $i = 0; $i < 5; ++$i){
	   my $aux_name = $names[0]."_".$i;
	   my $aux_seq = $alignment[0];
	   push @names, $aux_name;
 	   push @alignment, $aux_seq; 
     	}
    }
}

# Generates a file containing unaligned FASTA sequences for the validation step
my $raw_file = $output."/".$outfile_name."_validation_file.fasta";
open(FASTA, ">$raw_file");
for(my $i = 0; $i < scalar(@names); ++$i){
    my $aux = $alignment[$i];
    $aux =~ s/-//gi;
    $aux =~ s/\*//gi;
    print FASTA ">$names[$i]\n$aux\n";
}
close(FASTA);

# Removes gap-only columns from the aligment

my @matrix = @{$aux3};
my $size = scalar(@matrix);
my $num_cols = scalar(@{$matrix[0]});
my $count;
my $sub = 0;
my $c = 0;
for(my $i = 0; $i < $num_cols; ++$i){
    $count = 0;
    for(my $j = 0; $j < $size; ++$j){
    	if($matrix[$j][$i] eq '-'){
            ++$count;
        }
        else{
            last;
        }
     }
     if($count == $size){
        ++$c;
        my $pos = $i - $sub;
        for(my $j = 0; $j < $size; ++$j){
      	    my @aux = split("", $alignment[$j]);
            splice @aux, $pos, 1;
            $alignment[$j] = join("", @aux);
        }
        ++$sub;
    }
}

if($c == 0){
    print LOG "No gap-only sequences found.\n";
}
else{
    print LOG "$c gap-only columns removed.\n";
}

$no_red .= "_no_gap_only.fasta";
open(OUT, ">$no_red");
for(my $i = 0; $i < scalar(@names); ++$i){
    print OUT ">$names[$i]\n$alignment[$i]\n";
}

close(OUT);

my $seq_len = length($alignment[0]);
my $i;
my $window_score;
my $aux_len;

print LOG "\nPos-processing:\n";
print LOG "-clean Discard auxiliary files generated during the execution of Tabajara = $clean.\n";

for($i = 1; $i < scalar(@alignment); ++$i){
    $aux_len = length($alignment[$i]);
    if($aux_len != $seq_len){
	die "ERROR: Sequences of different lengths: $names[0] $aux_len != $names[$i] $seq_len\n";
    }
}

# Assigns weight for the alignment (if use_seq_weights option is true)

if(lc($use_seq_weights) eq 'true'){
    my $wei_file = $input_file;
    $wei_file =~ s/$align_suffix/weights/i;
    if(-e $wei_file){
	my $aux = load_sequence_weights($wei_file);
	@seq_weights = @{$aux};		
    }
    else{
	my $aux = calculate_sequence_weights(@alignment);
	@seq_weights = @{$aux};
    } 
}

if(scalar(@seq_weights) != scalar(@alignment)){
    @seq_weights = ();
    $len = scalar(@alignment);
    for($i = 1; $i < $len; ++$i){
	$seq_weights[$i] = 1 * $len;
    }
}

my @scores = ();
my @types = ();
my @values_sh = ();
my @values_mi = ();
my @valid_blocks_mi_sh = ();
my @valid_blocks_mi = ();
my @valid_blocks_sh = ();
my @final_blocks_sh = ();
my @final_blocks_mi = ();
$start = undef;
$end   = undef;
my $count_t = 0;
my $tp;
my @aux_t;
my @final = ();

# Generates full length profile HMMs (execution in full-length mode)

if(lc($full_length) eq "yes"){
    $hmmbuild =~ s/\"//g;
    my $type;
    if($file_type == 1){
        $type == "--dna";
    }
    else{
        $type = "--amino";
    }
    if(lc($measure) eq "c"){ #Conservation mode
        my $hmm_dir = $output."/hmm";
        if(!-e $hmm_dir){
            system "mkdir $hmm_dir";
        }
        my @amount = ();
        if(defined $cat){
            for(my $i = 0; $i < scalar(@categories); ++$i){
                my $count_cat = 0;
                foreach my $n (@names){
                    $n =~ s/>//g;
                    my @aux = split('_', $n);
                    if(uc($aux[0]) eq uc($categories[$i])){
                        ++$count_cat;
                    }
                }
                push @amount, $count_cat;
           }
        }
        my $hmm = $hmm_dir."/".$outfile_name.".hmm";
	# Runs hmmbuild program to build the profile HMM
        my $resp = system "hmmbuild $type $hmmbuild $hmm $input_file > /dev/null";
        if($resp > 0){
            print LOG "ERROR: Could not run hmmbuild: hmmbuild $hmmbuild $hmm $input_file\n";
            close(LOG);
            die "ERROR: Could not run hmmbuild: hmmbuild $hmmbuild $hmm $input_file\n";
        }

	# Validation of the generated model
        hmmValidation_conservation($hmm_dir, scalar(@alignment), $raw_file, $file_type, \@categories, \@amount);
	
	if($clean eq "yes"){
            system "rm -rf $hmm_dir/excluded_HMMs";
    	}
    }
    elsif(lc($measure) eq "d"){ # Discrimination mode
        my $hmm_dir = $output."/hmms";
        if(!-e $hmm_dir){
            system "mkdir $hmm_dir";
        }
        my $blocks = $output."/fastas";
        system "mkdir $blocks";

	# Generating the MSA for each category
        
	for(my $i = 0; $i < scalar(@categories); ++$i){
            my $category = $categories[$i];
            if(scalar(@categories) > 1){
                $hmm_dir = $output."/hmms/".$category;
                system "mkdir $hmm_dir";
            }
            my $count_cat = 0;
            my $raw_cat_file = $blocks."/".$category.".fasta";
            my $hmm = $hmm_dir."/".$category.".hmm";
            open(FASTA, ">$raw_cat_file");
            for(my $j = 0; $j < scalar(@names); ++$j){
            	my $name = $names[$j];
                $name =~ s/>//g;
                my @aux = split('_', $name);
                if(uc($aux[0]) eq uc($category)){
                   ++$count_cat;
                   my $seq = $alignment[$j];
                   $seq =~ s/-//gi;
                   $seq =~ s/\*//gi;
                   print FASTA ">$names[$j]\n$seq\n";
                }
            }
            close(FASTA);
	    if($count_cat > 0){
            	my $aln_cat_file = $blocks."/".$category."_aln.fasta";
            	my $resp = system "muscle -in $raw_cat_file -out $aln_cat_file 1> /dev/null 2> /dev/null";
            	if($resp > 0){
                    print LOG "ERROR: Could not run muscle: muscle -in $raw_cat_file -out $aln_cat_file\n";
                    close(LOG);
		    die "ERROR: ERROR: Could not run muscle: muscle -in $raw_cat_file -out $aln_cat_file\n";
            	}
            	my $resp = system "hmmbuild $type $hmmbuild $hmm $aln_cat_file > /dev/null";
            	if($resp > 0){
                    print LOG "ERROR: Could not run hmmbuild: hmmbuild $hmmbuild $hmm $aln_cat_file\n";
                    close(LOG);
                    die "ERROR: Could not run hmmbuild: hmmbuild $hmmbuild $hmm $aln_cat_file\n";
            	}

		# Validation of the generated models
            	hmmValidation($hmm_dir, $category, $count_cat, $raw_file, $file_type);
	    }
	    else{
		print LOG "\n$category sequences have not been found in the dataset. Unable to run in Discrimination Mode.\n";
		print STDERR "\n$category sequences have not been found in the dataset. Unable to run in Discrimination Mode.\n";
	    }
        }
    }
}
else{ 
    # Generates short profile HMMs

    if(lc($measure) eq "d"){ # Discrimination mode
    	my $blocks = $output."/blocks";
    	system "mkdir $blocks";
    	my $fastas_dir = $blocks."/original";
    	my $aux_original = $fastas_dir;
    	system "mkdir $fastas_dir";
    	my $without_red = $blocks."/final";
    	system "mkdir $without_red";
    	my $aux_final = $without_red;
    	my $hmms_dir = $output."/hmms";
    	system "mkdir $hmms_dir";
    	my $num_cat = scalar(@categories);
    	for(my $i = 0; $i < $num_cat; ++$i){
            $category = $categories[$i];
            if($num_cat > 1){
	    	print STDERR "\nCategory: $category\n";
            	$fastas_dir  = $blocks."/original/".$category;
            	system "mkdir $fastas_dir";
            	$without_red = $blocks."/final/".$category;
            	system "mkdir $without_red";
            	$hmms_dir    = $output."/hmms/".$category;
            	system "mkdir $hmms_dir";
            	print LOG "\n\nCategory: $category\n";
            }

    	    # Blocks selection by Mutual Information (MI) method

    	    print LOG "\nAnalysis by MI: \n";
    	    my $fastas_mi = $fastas_dir."/MI"; 
    	    system "mkdir $fastas_mi";	
	    my $aux_scores_mi;
	    my $aux_gaps_only;

	    # Calculates scores for all positions using MI method
    	    
	    ($aux_scores_mi, $aux_gaps_only) = mi(\@names, \@alignment, $gap_cutoff, $category);
	    if($aux_scores_mi == -1){
	    	next;
	    }
    	    my @scores_mi = @{$aux_scores_mi};
	    my @gaps_only = @{$aux_gaps_only};

	    # Selects the best blocks
	    my ($aux11, $aux12) = select_valid_blocks_disc(\@scores_mi, \@gaps_only, $threshold, $window_size, 1, $num_cat);

    	    @valid_blocks_mi = @{$aux12};
    	    @values_mi = @{$aux11};
    	    my $t = scalar(@valid_blocks_mi);
	    my $n = 0;
	    my $red_mi = $without_red."/MI";
            system "mkdir $red_mi";
	
	    # Analysis of the selected blocks

	    if($t > 0){
    	    	if(defined $interval){
            	    my $aux1 = interval_union(\@valid_blocks_mi, $interval, \@values_mi);
            	    @final_blocks_mi = @{$aux1};
    	    	}
    	    	else{
            	    @final_blocks_mi = @valid_blocks_mi;
    	    	}

            	@final_blocks_mi = @{eliminate_zeros(\@final_blocks_mi, \@scores_mi, $min_size_block)};
    	    	@valid_blocks_mi = ();
    	    	@valid_blocks_mi = @{createOriginalBlocksFiles($fastas_mi, \@final_blocks_mi, \@alignment, $measure, $gap_perc)};
		
    	    	if($max_block_size){
            	    createFinalBlocksSelectingBestRegion($fastas_mi, \@values_mi, $min_size_block, $max_block_size, $measure, $category, $input_file_prefix, $red_mi, $category, \@gaps_only);
    	    	}
    	    	else{
            	    createFinalBlocksByOriginal($fastas_mi, $red_mi, \@values_mi, $min_size_block);
    	    	}

	    	opendir(diretorio, "$red_mi");
    	    	my @fastas = readdir(diretorio);
    	    	closedir(diretorio);
   	    	foreach my $k (@fastas){
	    	    if($k eq "." or $k eq ".."){}
	    	    else{
		    	++$n;		
	    	    }
	    	}
	    }
            if($t == 0){
            	print STDERR "No block selected by MI!\n";
            	print LOG "No block selected by MI!\n\n";
            }
  	    else{
	    	print STDERR "Blocks selected by MI: $n \n";
            	print LOG "Blocks selected by MI: $n \n\n";
	    }	

	    # Blocks selection by Sequence Harmony (SH) method
	    
            print LOG "Analysis by SH: \n\n";
    	    my $fastas_sh = $fastas_dir."/SH";
	    system "mkdir $fastas_sh";
	    my $aux_scores_sh;
	    my $aux_gaps_only;

	    # Calculates scores for all position using SH method
    	    
	    ($aux_scores_sh, $aux_gaps_only) = sh(\@names, \@alignment, $gap_cutoff, $category);
	    my @gaps_only = @{$aux_gaps_only};
    	    my @scores_sh = @{$aux_scores_sh};

	    # Selecting the best blocks
    	    
	    my ($aux21, $aux22) = select_valid_blocks_disc(\@scores_sh, \@gaps_only, $threshold, $window_size, 2, $num_cat);
    	    @valid_blocks_sh = @{$aux22};
    	    @values_sh = @{$aux21};
    	    $t = scalar(@valid_blocks_sh);
	    my $red_sh = $without_red."/SH";
            system "mkdir $red_sh";

	    # Analysis of the selected blocks

	    if($t > 0){
            	if(defined $interval){
            	    my $aux1 = interval_union(\@valid_blocks_sh, $interval, \@values_sh);
            	    @final_blocks_sh = @{$aux1};
    	    	}
    	    	else{
            	    @final_blocks_sh = @valid_blocks_sh;
    	    	}

    	    	@final_blocks_sh = @{eliminate_zeros(\@final_blocks_sh, \@scores_sh, $min_size_block)};
 	    	@valid_blocks_sh = ();
    	    	@valid_blocks_sh = @{createOriginalBlocksFiles($fastas_sh, \@final_blocks_sh, \@alignment, $measure, $gap_perc)};
    	    	if($max_block_size){
            	    createFinalBlocksSelectingBestRegion($fastas_sh, \@values_sh, $min_size_block, $max_block_size, $measure, $category, $input_file_prefix, $red_sh, $category, \@gaps_only);
    	    	}
    	    	else{
            	    createFinalBlocksByOriginal($fastas_sh, $red_sh, \@values_sh, $min_size_block);
    	    	}

   	    	opendir(diretorio, "$red_sh");
            	my @fastas = readdir(diretorio);
            	closedir(diretorio);
            	$n = 0;
            	foreach my $k (@fastas){
            	    if($k eq "." or $k eq ".."){}
            	    else{
                    	++$n;
            	    }
            	}
	    }
	    if($t == 0){
            	print STDERR "No block selected by SH!\n";
            	print LOG "No block selected by SH!\n\n";
            }
            else{
            	$t = scalar(@final_blocks_sh);
            	print STDERR "Blocks selected by SH: $n\n";
            	print LOG "Blocks selected by SH: $n\n\n";
            }

    	    # Collapsing results obtained by MI and SH

	    @values = ();
	    @scores = ();
    	    for(my $i = 0; $i < scalar(@values_sh);++$i){
            	$values[$i] = ($values_mi[$i]+$values[$i])/2;
            	$scores[$i] = ($scores_mi[$i]+$scores_sh[$i])/2;
    	    }      

    	    my ($a1, $a2) = blocks_union($red_mi, $red_sh, $interval,\@values);
	    @valid_blocks = ();
    	    @valid_blocks = @{$a1}; 
   	    @types = @{$a2};     
	    my $t = scalar(@valid_blocks);
    	    print LOG "Generated region by union of MI and SH blocks:\n\n";
	    my $num_seq = 0;
    	    for(my $i = 0; $i < scalar(@valid_blocks); ++$i){ 	   
	    	++$num_seq;# = $i + 1;
            	my @aux = split("-", $valid_blocks[$i]);
            	my $t = "$aux[0]-$aux[1]";
	    	print LOG "Block: $t\n";
	    	if($types[$i] eq '1'){
	    	    my $file = `ls $red_mi/*_$t.fasta*`;
	    	    chomp($file);
     	    	    my $new_file = $without_red."/".$category."_".$num_seq."_".$t.".fasta";
	    	    system "cp $file $new_file";
	    	    if($max_block_size){
	   	    	print LOG "Top score region: $aux[0] - $aux[1]\n\n";
	    	    }
	    	}
	    	elsif($types[$i] == '2'){
            	    my $file = `ls $red_sh/*_$t.fasta*`;
            	    chomp($file);
	    	    my $new_file = $without_red."/".$category."_".$num_seq."_".$t.".fasta";
            	    system "cp $file $new_file";
	    	    if($max_block_size){
                    	print LOG "Top score region: $aux[0] - $aux[1]\n\n";
            	    }
	    	}
	    	elsif($types[$i] == 3){
	    	    my $new_file = $without_red."/".$category."_".$num_seq."_".$t.".fasta";
	    	    my @new_names = ();
	    	    my @new_seqs = ();
	    	    for(my $j = 0; $j < scalar(@names); ++$j){
		    	if($names[$j] =~ /$category/){
		    	    my $len = ($aux[1]-1)-($aux[0]-1)+1;
		    	    push @new_names, $names[$j];
		    	    my $seq = substr($alignment[$j], ($aux[0]-1), $len);
		    	    push @new_seqs, $seq;
		    	}
	    	    }
	    	    my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@new_names, \@new_seqs, undef);
            	    @new_names = @{$aux_name};
            	    @new_seqs = @{$aux_seqs};
	    	    my $min_nseq = 2;
	    	    my $t = scalar(@new_names);
	    	    if(defined $amount){
                    	$min_nseq = $amount;
            	    }   
            	    if(scalar(@new_names) < $min_nseq){
		    	print LOG "Number of sequences lower than $min_nseq: discarded\n\n";
                    	next;
            	    }
	    	    my $len = length($new_seqs[0]);
	    	    my $size = scalar(@new_seqs);
	    	    my @delete = ();
	    	    for(my $k = 0; $k < $len; ++$k){
	    	    	my @col = @{get_column($k, \@new_seqs)};
	    	    	my $gap = 100*(gap_percentage(\@col));
	    	    	if($gap == 100){
		    	    for(my $l = 0; $l < $size; ++$l){
		            	push @delete, ($aux[0]-1+$k); 
		    	    }
		    	    --$len;
	    	    	}
	    	    }
	    	    $len = length($new_seqs[0]);
	    	    if($len < $min_size_block){
	    	    	print LOG "Selected region is shorter than the minimum block size after removing gap-only columns: discarded\n";
	    	    	next;
	    	    }
	    	    if($max_block_size){
		    	my $size = $aux[1] - $aux[0] + 1;
        	    	if($size >= $max_block_size){
		    	    my @blocks = ();
            	    	    for(my $i = 0; $i < $size; ++$i){
                	    	my $count = 0;
                	    	my $len_del = scalar(@delete);
                	    	if($len_del > 0){
                    	    	    for(my $k = 0; $k > $len_del; ++$k){
                        	    	if($delete[$k] >= $aux[0] and $delete[$k] <= $aux[1]){
                            	    	    $count++;
                        	    	}
                            	    }
                	    	}
                	    	my $s = $i;
                	    	my $e = ($i + $max_block_size+$count)-1;
                	    	if($e >= $size){
                    	    	    last;
                	    	}
                	    	my $sum = 0;
                	    	my $aux_start = $aux[0] - 1 + $i;
                	    	my $aux_end = $aux[0] - 1  + $e;
			    	my $win = 0;
                	    	for(my $j = $aux_start; $j < $aux_end; ++$j){
				    if($win >= $max_block_size){
                		    	last;
            			    }
				    if($aux_end > $aux[1]){
				    	last;
				    }
            			    if($gaps_only[$j] == 1){
                		    	++$aux_end;
            			    }
            			    elsif($values[$j] > 0){ #!= -1000){
                		    	$sum += $values[$j];
                		    	++$win;
            			    }
                	    	}
                	    	my $string = "$sum-$s-$e";
                	    	my $as = $aux_start + 1;
                	    	my $ae = $aux_end + 1;
                	    	print LOG "Region: $as - $ae: score sum=$sum\n";
                	    	push @blocks, $string;
		    	    }
		    	    my @aux2 = split("-", $blocks[0]);
            	    	    my $sum_max = $aux2[0];
            	    	    my $start_max = $aux2[1];
            	    	    my $end_max = $aux2[2];

            	    	    for(my $k = 1; $k < scalar(@blocks); ++$k){
                	    	my @aux2 = split("-", $blocks[$k]);
                	    	if($aux2[0] > $sum_max){
                    	    	    $sum_max = $aux2[0];
                    	    	    $start_max = $aux2[1];
                    	    	    $end_max = $aux2[2];
                	    	}
            	    	    }
		    	    my $as = $start_max + $aux[0];
            	    	    my $ae = $end_max + $aux[0];
            	    	    my ($a, $e) = trimming_zeros(($as-1), ($ae-1), \@values, $min_size_block);
            	    	    print LOG "Top score region: $as - $ae: score sum=$sum_max\n";
		    	    if(($e - $a + 1) >= $min_size_block){
                	    	my $aux_start = $a - $start;
                	   	my $aux_end = $e - $start;
                	    	++$a;
                	    	++$e;
			    	print LOG "Selected region after zero-score trimming:$a-$e\n";
			    	my @new_seqs2 = ();
			    	my $start_a = ($a - 1) - ($aux[0] - 1);
                            	my $end_a = ($aux[1] - 1) - ($e - 1);
                            	my $len = ($e - 1) - ($a - 1) + 1;
                	    	for(my $k = 0; $k < scalar(@new_names); ++$k){
                    	    	    my $str = substr $new_seqs[$k], $start_a, $len;
                    	    	    $new_seqs2[$k] = $str;
                	    	}
			        my $len = length($new_seqs2[0]);
                	        my $size = scalar(@new_seqs2);
                	        for(my $k = 0; $k < $len; ++$k){
                    	    	    my @col = @{get_column($k, \@new_seqs2)};
                    	    	    my $gap = 100*(gap_percentage(\@col));
                    	    	    if($gap == 100){				
                        	    	for(my $l = 0; $l < $size; ++$l){				    
                            	    	    my @aux = split("", $new_seqs2[$l]);
                            	    	    splice @aux,$k, 1;
                            	    	    $new_seqs2[$l] = join("", @aux);
                        	    	}
                        	    	--$len;
                    	    	    }
                	    	}
                	    	$len = length($new_seqs2[0]);
                	    	if($len < $min_size_block){
                    	    	    print LOG "Selected region is shorter than the minimum block size after removing gap-only columns: discarded\n";
                    	    	    next;
                	    	}
			    	my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@new_names, \@new_seqs2, undef);
                	    	@new_names = @{$aux_name};
                	    	@new_seqs2 = @{$aux_seqs};
			    	my $nseq = scalar(@new_names); 
        		    	my $min_nseq = 2;
        		    	if(defined $amount){
            	 	    	    $min_nseq = $amount;
        		    	}
        		    	if($nseq < $min_nseq){
            		    	    print LOG "Number of sequences lower than $min_nseq: discarded\n\n";
            		    	    next;
        		    	}

			    	$new_file = $without_red."/".$category."_".$num_seq."_".$as."-".$ae.".fasta";
                	    	open(FILE, ">$new_file");
                	    	for(my $k = 0; $k < scalar(@new_names); ++$k){
                    	    	    print FILE ">$new_names[$k]\n";
                    	    	    print FILE "$new_seqs2[$k]\n";
                	    	}
                	    	close(FILE);
		    	    }
		    	    
			    print LOG "\n";
			    
			    # Tries to select other regions 
                	    
			    if(lc($more_blocks) eq "yes"){
                    	    	print LOG "Trying to select more subregions:\n";
                    	    	my @aux_blocks = @{searchForMoreRegions($aux[0], $aux[1], ($as - 1), ($ae - 1), \@values, \@gaps_only)};
                    	    	if(scalar(@aux_blocks) == 0){
                        	    print LOG "No new subregion found!\n\n";
                    	    	}
                    	    	else{
                        	    my $last = $ae - 1;
                        	    my $first = $as - 1;
				    my @names_cat = ();
                        	    my @aln_cat = ();
                        	    for(my $c = 0; $c < scalar(@names); ++$c){
                           	    	my $name = $names[$c];
                            	    	$name =~ s/>//g;
                            	    	my @aux = split('_', $name);
                             	    	if(uc($aux[0]) eq uc($category)){
                                	    push @names_cat, $names[$c];
                                	    push @aln_cat, $alignment[$c];
                            	    	}
                        	    }
                        	    @aux_blocks = @{sortCoordinates(\@aux_blocks)};
                        	    for(my $l = 0; $l < scalar(@aux_blocks); ++$l){
                            	    	my @aux_b = split("-", $aux_blocks[$l]);
                            	    	my $start_b = $aux_b[0];
                            	    	if(($start_b == $last) or ($start_b == $e - 1)){
                                	    ++$start_b;
                            	        }
                            	    	elsif(($start_b < $last) and (($last - $start_b) < 5)){
                                	    my $dif = $last-$start_b;
                                	    $start_b += ($dif + 1);
                            	    	}
                            	    	my $end_b = $aux_b[1];
                            	    	if(($end_b == ($a - 1)) or ($end_b == $first)){
                                	    --$end_b;
                            	    	}
                            	    	if(($end_b - $start_b + 1) < $min_size_block){
                                	    next;
                            	    	}
                            	    	my $aux_s = $start_b + 1;
                            	    	my $aux_e = $end_b + 1;
				    	my ($ts, $te) = trimming_gaps_only($start_b, $end_b, \@gaps_only, $min_size_block);
                            	    	print LOG "New block: $ts - $te\n";
                            	    	++$num_seq;
                            	    	my $resp = generateRegion($start, $end, $ts, $te, \@values, \@names_cat, \@aln_cat, $without_red, $num_seq, $category);
                            	    	if($resp == -1){
                                	    --$num_seq;
                            	    	}
                            	    	$last = $end_b;
                            	    	$first = $start_b;
                        	    }
		    	    	}
		    	    }
		    	}
		    	elsif($size < $min_size_block){
		    	    print LOG "Selected region is shorter than the minimum block size: discarded\n\n";
		    	}
		    	else{
		    	    my $len = length($new_seqs[0]);
                    	    my $size = scalar(@new_seqs);
                    	    for(my $k = 0; $k < $len; ++$k){
                    	    	my @col = @{get_column($k, \@new_seqs)};
                            	my $gap = 100*(gap_percentage(\@col));
                            	if($gap == 100){
                            	    for(my $l = 0; $l < $size; ++$l){
                            	    	my @aux = split("", $new_seqs[$l]);
                                    	splice @aux,$k, 1;
                                    	$new_seqs[$l] = join("", @aux);
                            	    }
                            	    --$len;
                            	}
                      	    }
                    	    $len = length($new_seqs[0]);
                    	    if($len < $min_size_block){
                            	print LOG "Selected region is shorter than the minimum block size after removing gap-only columns: discarded\n";
                            	next;
                    	    }
                    	    my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@new_names, \@new_seqs, undef);
                    	    @new_names = @{$aux_name};
                    	    @new_seqs = @{$aux_seqs};
		    	    my $nseq = scalar(@new_names); 
        	    	    my $min_nseq = 2;
        	    	    if(defined $amount){
            		    	$min_nseq = $amount;
        	    	    }
        	    	    if($nseq < $min_nseq){
            		    	print LOG "Number of sequences lower than $min_nseq: discarded\n\n";
            		    	next;
        	    	    }

                    	    open(FILE, ">$new_file");
                    	    for(my $k = 0; $k < scalar(@new_names); ++$k){
                    	    	print FILE ">$new_names[$k]\n";
                            	print FILE "$new_seqs[$k]\n";
                      	    }
                    	    close(FILE);
		    	}
	    	    }
	    	    else{
	    	    	open(FILE, ">$new_file");
	    	    	for(my $k = 0; $k < scalar(@new_names); ++$k){
                    	    print FILE ">$new_names[$k]\n";
            	    	    print FILE "$new_seqs[$k]\n";
            	    	} 
	    	    	close(FILE);
	    	   }
	    	}
    	    }
    	    system "rm -rf $red_mi $red_sh";
	    my $qnt = `ls $without_red | wc -l`; 
	    print STDERR "Blocks selected by MI+SH: $qnt\n";
            print LOG "Blocks selected by MI+SH: $qnt\n";
	    
	    # Profile HMM construction
	    my ($aux_final, $msg) = createHMMs($without_red, $hmms_dir, $measure, $category);
	    @final = @{$aux_final};
	    if(defined $msg){
	    	print LOG $msg;
	    	print STDERR $msg;
	    	exit;
	    }
	    my $count_cat = 0;
            foreach my $n (@names){
            	$n =~ s/>//g;
            	my @aux = split('_', $n);
            	if(uc($aux[0]) eq uc($category)){
                    ++$count_cat;
            	}
            }

	    # Profile HMM validation
	    
	    @final = @{hmmValidation($hmms_dir, $category, $count_cat, $raw_file, $file_type)};
    	    @final = sort { $a <=> $b } @final;
	    $t = scalar(@final);
	    print STDERR "Total # of valid blocks: $t\n";
	    print LOG "Total # of valid blocks: $t\n";
	    if($num_cat > 1){
    	    	generateFinalFilesMix(\@final, \@final_blocks_mi, \@final_blocks_sh, \@values, \@values_mi, \@values_sh, \@scores, $blocks, $category);
	    }
	    else{
	    	generateFinalFilesMix(\@final, \@final_blocks_mi, \@final_blocks_sh, \@values, \@values_mi, \@values_sh, \@scores, $blocks, undef);
	    }   
    	}
    	if($clean eq "yes"){
            system "rm -rf $aux_original $aux_final";
    	}
    }
    elsif(lc($measure) eq "c"){ # Conservation mode
    	my $blocks = $output."/blocks";
    	system "mkdir $blocks";
    	my $fastas_dir = $blocks."/original"; 
    	my @gaps_only = ();
    	if($original_size == 1){
	    @scores = ();
	    my $len = length($alignment[0]);
	    for(my $i = 0; $i < $len; ++$i){
	    	$scores[$i] = 1;
	    }
    	}
    	else{
    	   if($file_type == 1){ # Nucleotide sequences

		# Calculates scores for all positions of the MSA using Shannon Entropy
	    	
		my $aux_scores = shannon_conservation(\@{$aux3}, $gap_cutoff);
    	    	@scores = @{$aux_scores};
    	   }
    	   else{ # Amino acid sequences
		
		# Calculates scores for all positions of the MSA using Jensen-Shannon divergency
    	   	
		@scores = ();   
	    	my @matrix = @{$aux3};
	    	my $size = scalar(@matrix);
	    	for($i = 0; $i < (length($alignment[0])); ++$i){
		    my @col_matrix = ();
		    my $gap = 0;
		    for(my $j = 0; $j < $size; ++$j){
                    	$col_matrix[$j] = $matrix[$j][$i];
		    	if($matrix[$j][$i] eq "-"){
                            ++$gap;
                    	}
                    }
		
		    my $t = scalar(@col_matrix);
		    my $p = $size;
                    if(scalar(@col_matrix) == scalar(@alignment)){
		    	$gap = ($gap/$size)*100;
                    	if($gap <= $gap_cutoff){
                            my $score_position = js_divergence(\@col_matrix, \@bg_distribution, \@seq_weights, $use_gap_penalty );
                            push @scores, $score_position;
                    	}
                   	else{
                            push @scores, (-1000);
                    	}
                    }
            	}

    	    	# Updates scores for window size
    	    	
		if($window_size_score > 0){
    	    	    my $aux_scores = window_score(\@scores, $window_size_score, $win_lam);
    	    	    @scores = ();
    	    	    @scores = @{$aux_scores};
    	    	}  
    	   }
    	}
        # Normalizes scores if the normalize_scores option is true
        
	# Selection of the best blocks
        
	my ($aux1, $aux2) = select_valid_blocks(\@scores, $threshold, $window_size, 3, 1);
        @valid_blocks = @{$aux2};
        @values = @{$aux1};
    
        if(scalar(@valid_blocks) == 0){
	    print LOG "No selected blocks!\n";
	    print STDERR "No selected blocks!\n";
  	    my $file = $blocks."/scores.csv";
            open(FILE, ">$file");
            print FILE "Position\tAll scores\n";
            for($i = 0; $i < scalar(@values); ++$i){
            	my $p = $i+1;
            	if($values[$i] == -1000){
                    print FILE "$p\t0\n";
            	}
            	else{
                    print FILE "$p\t$values[$i]\n";
            	}
            }
            close(FILE);
            close(LOG);
	    if(lc($clean) eq "yes"){
	    	my $num = `ls $output/*.fasta | wc -l`;
	    	if($num > 0){
		    system "rm $output/*.fasta";
	    	}
   	    } 	
	    exit;
    	}
    	if(defined $interval){
    	    my ($aux1) = interval_union(\@valid_blocks, $interval,\@values);
            @final_blocks = @{$aux1};
    	}
    	else{
            @final_blocks = @valid_blocks;
    	}
    	@final_blocks = @{eliminate_zeros(\@final_blocks, \@scores, $min_size_block)};
    	$t = scalar(@final_blocks);

    	@valid_blocks = ();
    	system "mkdir $fastas_dir";
    	@valid_blocks = @{createOriginalBlocksFiles($fastas_dir, \@final_blocks, \@alignment, $measure, $gap_perc)};
    	my $without_red = $blocks."/final";
    	system "mkdir $without_red";
    	if($max_block_size){
	    createFinalBlocksSelectingBestRegion($fastas_dir, \@values, $min_size_block, $max_block_size, $measure, $category, $input_file_prefix, $without_red,$input_file_prefix, \@gaps_only);
    	}
    	else{
            createFinalBlocksByOriginal($fastas_dir, $without_red, \@values, $min_size_block);
    	}
    	my @amount = ();
    	if(defined $cat){
	    for(my $i = 0; $i < scalar(@categories); ++$i){
	        my $count_cat = 0;
	        foreach my $n (@names){
    		    $n =~ s/>//g;
            	    my @aux = split('_', $n);
            	    if(uc($aux[0]) eq uc($categories[$i])){
                    	++$count_cat;
            	    }
	    	}
	    	$amount[$i] = $count_cat;
            }
    	}

        my $hmms_dir = $output."/hmms";
        system "mkdir $hmms_dir";

	# Profile HMM construction
	
        my ($aux_final, $msg) = createHMMs($without_red, $hmms_dir, $measure, $category);
        if(defined $msg){
            print LOG $msg;
            print STDERR $msg;
    	    exit;
    	}

	# Profile HMM validation
	 
    	@final = @{hmmValidation_conservation($hmms_dir, scalar(@alignment), $raw_file, $file_type, \@categories, \@amount)};
    	@final = sort { $a <=> $b } @final;
    	$t = scalar(@final);
    	print LOG "Blocks selected: $t \n\n";
    	print STDERR "Blocks selected: $t \n";

	# Generating the final files
	
    	generateFinalFiles(\@final, \@values, \@scores, $blocks, undef);
    	my $all_blocks = $output."/all_selected_blocks.fasta";
    	open(ALL, ">$all_blocks");
    	for(my $k = 0; $k < scalar(@names); ++$k){
	    print ALL ">$names[$k]\n";
	    my $str = "";
	    my $aux_seq = $alignment[$k];
	    for(my $b = 0; $b < scalar(@final); ++$b){
	    	my @aux_c = split("-", $final[$b]);
	    	my $start = $aux_c[0]-1;
	    	my $end = $aux_c[1]-1;
	    	my $len = $end-$start+1;
	    	$str .= substr $aux_seq, $start, $len;
	    }
	    print ALL "$str\n";
    	}
    	close(ALL);
    	if($clean eq "yes"){
            system "rm -rf $fastas_dir $without_red $hmms_dir/excluded_HMMs";
   	}
    }
} 

system "rm $raw_file";

if($clean eq "yes"){
    system "rm $no_red";
}

close(LOG);

################################################################################
#                                   Subroutines                                #
################################################################################ 

################################################################################
#                            Frequency Count and Gap Penalty                   #
################################################################################

sub weighted_freq_count_pseudocount{
    my ($col_aux, $wei_aux, $pc_amount) = @_;
    my @col = @{$col_aux};
    my @seq_weights = @{$wei_aux} ;
    my $i;
    my $j;
    
    if (scalar(@seq_weights) != scalar(@col)){
	my $len = scalar(@col);
	@seq_weights = ();
	for($i = 0; $i < $len; ++$i){
            $seq_weights[$i] = 1;
	}
    }
    my  $aa_num = 0;
    my @freq_counts = ();
    my $len_amino;
    if($file_type == 1){
    	$len_amino = scalar(@nucleotides);
    }
    else{
	$len_amino = scalar(@amino_acids);
    }
    for($i = 0; $i < $len_amino; ++$i){
	$freq_counts[$i] = $pc_amount;
    }
    for($i = 0; $i < $len_amino; ++$i){
	for($j = 0; $j < scalar(@col); ++$j){
	    if($amino_acids[$i] eq $col[$j]){
		$freq_counts[$aa_num] += 1 * $seq_weights[$j];
	    }
	}
	++$aa_num;
    }  
    my $sum_wei = sum_arrays(\@seq_weights);
    for($i = 0; $i < scalar(@freq_counts); ++$i){
	$freq_counts[$i] = $freq_counts[$i]/($sum_wei + $len_amino*$pc_amount);	
    }
    return \@freq_counts;

}


################################################################################
#                              Reads Clustal alignment                          #
################################################################################
sub read_clustal_alignment{
    my $file = shift;
    my @nms = ();
    my @aligns = ();
    my %hash;
    my $count = 0;
    my %aux;
    my @matrix = ();
    open(FILE, "$file") or die "Could not open the $file file\n";
    while(<FILE>){
        chomp($_);
	if($_ =~ /\*/ or $_ =~ /\./){
	    next;
	}
	else{
	    my @aux = split(" ", $_);
	    if(scalar(@aux) == 2){
		my $aln = $aux[1];
		$aln = uc $aln;
		$aln =~ s/\r//gi;
		if(defined $hash{$aux[0]}){
		    $hash{$aux[0]} .= $aln;
		}
		else{
		    $hash{$aux[0]} = $aln;
		    $aux{$aux[0]} = $count;
		    ++$count;
		}
	    }
	}
    }
    close(FILE); 
    my $id;
    my $count = 0;
    foreach my $key (keys %hash){
        $id = $aux{$key};
	$nms[$id] = $key;
	$aligns[$id] = $hash{$key};
	my @aux = split("", $aligns[$id]);
	for(my $j = 0; $j < scalar(@aux); ++$j){
            $matrix[$count][$j] = $aux[$j];
        }
	++$count;
    }

    return (\@nms, \@aligns, \@matrix);
}

################################################################################
#                               Read FASTA alignment                           #
################################################################################
sub read_fasta_alignment_short{
    my $file = shift;
    my @nms = ();
    my @aligns = ();
    my $seq = undef;
    open(FILE, "$file") or die "Could not open the $file file\n";
    while(<FILE>){
        chomp($_);
	if($_ =~ />/){
	    my $name = $_;
	    $name =~s/>//i;
	    push @nms, $name;
	    if(defined $seq){
		$seq = uc $seq;
                $seq =~ s/\r//gi;
		push @aligns, $seq;
	    }
	    $seq = "";
	}
	else{
	    $seq .= $_;
	}
    }
    $seq = uc $seq;
    $seq =~ s/\r//gi;
    push @aligns, $seq;
    close(FILE);
    return (\@nms, \@aligns);
}

sub read_fasta_alignment{
    my $file = shift;
    my @nms = ();
    my @aligns = ();
    my $seq = undef;
    my @matrix = ();
    my $count = 0;
    open(FILE, "$file") or die "Could not open the $file file\n";
    while(<FILE>){
        chomp($_);
        if($_ =~ />/){
            my $name = $_;
            $name =~s/>//i;
            push @nms, $name;
            if(defined $seq){
                $seq = uc $seq;
                $seq =~ s/\r//gi;
                my @aux = split("", $seq);
                for(my $j = 0; $j < scalar(@aux); ++$j){
                    $matrix[$count][$j] = $aux[$j];
                }
                push @aligns, $seq;
                ++$count;
            }
            $seq = "";
        }
        else{
            $seq .= $_;
        }
    }
    $seq = uc $seq;
    $seq =~ s/\r//gi;
    my @aux = split("", $seq);
    for(my $j = 0; $j < scalar(@aux); ++$j){
        $matrix[$count][$j] = $aux[$j];
    }
    push @aligns, $seq;
    close(FILE);
        return (\@nms, \@aligns, \@matrix);
}

################################################################################
#                             Loads sequence weights                           #
################################################################################
sub load_sequence_weights{
    my $file = shift;
    open(FILE, "$file") or die "Could not open the $file file\n";
    while(<FILE>){
        chomp($_);
    }
}

################################################################################
#                             Calculates sequence weights                      #
################################################################################
sub calculate_sequence_weights{
    my @msa = @_;
    my $length = scalar(@msa);
    my @weights = (); # * $length;
    my $i;
    my $j;
    my $k;
    for($i = 0; $i < $length; ++$i){
	$weights[$i] = 0;
    }
    my $msa_length = length(@msa[0]);
    my @freq_counts;
    my @aux = ();

    # Constructs a matrix from the alignment sequences
    
    for($i = 0; $i < $length; ++$i){
	my @seq = split("", $msa[$i]);
	for($j = 0; $j < $msa_length; ++$j){
	    $aux[$i][$j] = $seq[$j];
	}
    }
    for($j = 0; $j < $msa_length; ++$j){
	@freq_counts = ();
	for($k = 0; $k < scalar(@amino_acids); ++$k){
                $freq_counts[$k] = 0;
        }
	for($i = 0; $i < $length; ++$i){	
	    if(!($aux[$i][$j] eq '-')){
		 $freq_counts[$aa_to_index{$aux[$i][$j]}] += 1;
	    }		
	}
	my $num_observed_types = 0;
	for($i = 0; $i < scalar(@freq_counts); ++$i){
            if ($freq_counts[$i] > 0){
		$num_observed_types +=1;
	    }
	}
        for($i = 0; $i < $length; ++$i){
            my $d = $freq_counts[$aa_to_index{$aux[$i][$j]}] * $num_observed_types;
            if ($d > 0){
            	$weights[$i] += 1/$d;
	    }
	}
    }
    for($i = 0; $i < scalar(@weights); ++$i){
    	$weights[$i] /= $msa_length;
    }
    return \@weights;
}

################################################################################
#             Returns the index of an element contained in an array            #
################################################################################
sub index_of{
    my ($array, $element) = @_;
    my @aux = @{$array};
    my $i;
    my $index = -1;
    for($i = 0; $i < scalar(@aux); ++$i){
  	if($aux[$i] eq $element){
	   $index = $i;
	}
    }
    return $index;
}

################################################################################
#                Returns the elements of an alignment column                   #
################################################################################
sub get_column{
    my ($index, $array) = @_;
    my @aux = @{$array};
    my $i;
    my @elements = ();
    for($i = 0; $i < scalar(@aux); ++$i){
	my @seq = split("", $aux[$i]);
	if($index < scalar(@seq)){
	    push @elements, $seq[$index];
	}
    }
    return \@elements;
}

################################################################################
#                Calculates the gap percentage of a column                     #
################################################################################
sub gap_percentage{
    my $aux = shift;
    my @col = @{$aux};
    my $gap = 0;
    my $i;
    for($i = 0; $i < scalar(@col); ++$i){
	if($col[$i] eq '-'){
	    ++$gap;
	}
    }
    return ($gap/scalar(@col));
}

################################################################################
#           Calculates of Jensen-Shannon Divergence for each column            #
################################################################################
sub js_divergence{
    my ($aux_col, $aux_dist, $aux_wei, $gap_penalty) = @_;
    my @distr = @{$aux_dist};
    my @seq_weights = @{$aux_wei};
    my @col = @{$aux_col};
    my $aux_fc = weighted_freq_count_pseudocount(\@col, \@seq_weights, PSEUDOCOUNT);
    my @fc = @{$aux_fc};
    my $i;

    if(scalar(@distr) == 20){
	my @new_fc = @fc;
	pop @new_fc;
	my $sum = sum_arrays(\@new_fc);

	for($i = 0; $i < scalar(@new_fc); ++$i){
	   $new_fc[$i] = $new_fc[$i]/$sum;
	}
	@fc = ();
	@fc = @new_fc;
    }
    my $t = scalar(@fc);
    my $k = scalar(@distr);
    if(scalar(@fc) != scalar(@distr)){
	return -1;
    }

    my @r = ();
    for($i = 0; $i < scalar(@fc); ++$i){
	$r[$i] = (0.5*$fc[$i])+(0.5*@distr[$i]); 
    }

    my $d = 0;

    for($i = 0; $i < scalar(@fc); ++$i){
	if($r[$i] != 0){
	    if($fc[$i] == 0){
		$d += $distr[$i] * log2($distr[$i]/$r[$i]);
	    }
	    elsif($distr[$i] == 0){
		$d += $fc[$i] * log2($fc[$i]/$r[$i]);
	    }
	    else{
		$d += ($fc[$i] * log2($fc[$i]/$r[$i])) + ($distr[$i] * log2($distr[$i]/$r[$i]));
	    }
	}
    }

    $d = $d/2;
  
    if($gap_penalty == 1){
	return ($d * weighted_gap_penalty(\@col, \@seq_weights));
    }
    else{
   	return $d;
    }
    
}

################################################################################
#        Calculates the Shannon entropy for each column of the alignment       #
################################################################################
sub shannon_conservation{
   my $aux = shift;
   my $gap_cut = shift;
   my @entropy = ();
   my @matrix = @{$aux};
   my $size = scalar(@matrix);
   my $num_cols = scalar(@{$matrix[0]});
   for(my $i = 0; $i < $num_cols; ++$i){
   	my %hash = %{initialize_DNA()};
        for(my $j = 0; $j < $size; ++$j){
             if(defined $hash{$matrix[$j][$i]}){
             	++$hash{$matrix[$j][$i]};
             }
        }
        my $ent = 0;	
  	if(($hash{"-"}/$size)*100 > $gap_cut){
	     push @entropy, -1000;
	}
	else{
             foreach my $key (sort keys %hash){
		if(lc($key) eq 'x'){}
		else{ 
             	    my $pos = ($hash{$key}/$size);
             	    if($pos > 0){
             	    	$ent += (($pos*log2($pos))/log2(5));
             	    }
		}
             }
             $ent *= -1;
             push @entropy, $ent;
	}
    }
    my $max_ent = -5*((1/5)*log2(1/5));
    for(my $i = 0; $i < scalar(@entropy); ++$i){
        if($entropy[$i] == -1000){
        }
        else{
            $entropy[$i] = 1 - $entropy[$i];
        }
    }
    return \@entropy;
}

################################################################################
#                           Calculates log2 of a number                        #  
################################################################################
sub log2{
    my $number = shift;
    my $result = log($number)/log(2);

    return $result;
}

################################################################################
#                           Calculates weighted gap penalty                    #
################################################################################
sub weighted_gap_penalty{
    my ($aux_col, $aux_wei) = @_;
    my @col = @{$aux_col};
    my @seq_weights = @{$aux_wei};
    my $i;	
    my $len_col = scalar(@col);
    if(scalar(@seq_weights) != $len_col){
	@seq_weights = ();
	for($i = 0; $i < $len_col; ++$i){
	    $seq_weigths[$i] = 1;
	}
    }
    my $gap_sum = 0;
    for($i = 0; $i < $len_col; ++$i){
        if ($col[$i] eq '-'){
            $gap_sum += $seq_weights[$i];
	}
    }
    my $sum = sum_arrays(\@seq_weights);
    my $result;
    if($sum == 0){
   	$result = 1;
    }
    else{
    	$result = 1 - ($gap_sum/sum_arrays(\@seq_weights)); 
    }
    return $result;

}

################################################################################
#                     Calculates the sum of array elements                     #
################################################################################
sub sum_arrays{
    my $aux_array = shift;
    my @array = @{$aux_array};
    my $i;
    my $sum = 0;
    for($i = 0; $i < scalar(@array); ++$i){
	$sum += $array[$i];
    }

    return $sum;
}

################################################################################
#                                Calculates window scores                      #
################################################################################
sub window_score{
    my ($aux_score, $window_len, $lam) = @_;
    my @w_scores = @{$aux_score};
    my $i;
    my $j;
    for($i = $window_len; $i <= (scalar(@w_scores)-$window_len); ++$i){
	if($scores[$i] < 0){
	    next;
	}
	my $sum = 0;
	my $num_terms = 0;
	for($j = $i - $window_len; $j <= $i + $window_len + 1; ++$j){
	    if($i != $j and $scores[$j] >= 0){
		$num_terms += 1;
		$sum += $w_scores[$j];
	    }
	}
	if($num_terms > 0){
	    $w_scores[$i] = (1 - $lam) * ($sum / $num_terms) + $lam * $w_scores[$i];
	}
    }
    return \@w_scores;
}

################################################################################
#                             Calculates z-scores                              #
################################################################################
sub calc_z_scores{
    #score_cutoff are not included.
    my ($aux_score, $cutoff) = @_;
    my @scores = @{$aux_score};

    my $average = 0;
    my $std_dev = 0;
    my @z_scores = ();
    my $num_scores = 0;
    my $i;
    
    for($i = 0; $i < scalar(@scores); ++$i){
	if($scores[$i] > $cutoff){
	    $average += $scores[$i];
	    $num_scores++;
	}
    }
    if($num_scores != 0){
	$average = $average/$num_scores;
    }
    
    for($i = 0; $i < scalar(@scores); ++$i){
	if($scores[$i] > $cutoff){
	   $std_dev += (($scores[$i] - $average)**2) / $num_scores;
	}
    }
    $std_dev = sqrt($std_dev);

    for($i = 0; $i < scalar(@scores); ++$i){
        if($scores[$i] > $cutoff and $std_dev != 0){
	    push @z_scores, (($scores[$i]-$average)/$std_dev);
	}
   	else{
	    push @z_scores, (-1000);
	}
    }
    return \@z_scores;
}

###################################################################################
# Checks if a directory with the given output name already exists. If so, a numeric
# suffix is added to the name.
###################################################################################

sub output_dir_name {
    my $output_dir_name = shift;
    my $count = 2;
    my $flag = 0;
    if (-d $output_dir_name) {
        $flag = 1;
        while (-d "$output_dir_name\_$count") {
            $count++;
        }
        $output_dir_name = "$output_dir_name\_$count";
    }
    print "\nOutput directory already exists, saving results to $output_dir_name instead.\n\n" if $flag;
    return ($output_dir_name);
}

####################################################################################
#                          Verifies the type of input file                         #
####################################################################################

sub verify_file_type{
    my $file = shift;
    open(FILE, $file);
    while(<FILE>){
	chomp($_);
	if($_ =~ />/){
	    close(FILE);
	    return 1;
	}
    }
    close(FILE);
    return 0;
}

######################################################################################
#         Calculates Mutual information for each position of the alignment           #
######################################################################################

sub mi{
    my $aux1 = shift;
    my $aux2 = shift; 
    my $gap  = shift; 
    my $category = shift;
   
    my @aux_n = @{$aux1};
    my @aux_a = @{$aux2};
    my $t = scalar(@aux_n);
    my @cat1;
    my @seq1;
    my @cat2;
    my @seq2;

    for(my $i = 0; $i < scalar(@aux_n); ++$i){
        my $name = $aux_n[$i];
        $name =~ s/>//g;
        my @aux = split('_', $name);
        if(uc($aux[0]) eq uc($category)){
            push @cat1, $name;
            push @seq1, $aux_a[$i];
        } 
        else{
            push @cat2, $name;
            push @seq2, $aux_a[$i];
        }
    } 
    my @names = (@cat1, @cat2);
    my @gaps_only = @{get_gapOnly(\@seq1)};

    my @alignment = (@seq1, @seq2);

    my $n1 = scalar(@cat1);
    my $n2 = scalar(@cat2);
    my $teste = scalar(@aux_n);	

    my $n_tot = $n1 + $n2;
    if($n1 == 0){
	print LOG "\n$category sequences have not been found in the dataset. Unable to run in Discrimination Mode.\n"; #Aborting...\n";
	print STDERR "\n$category sequences have not been found in the dataset. Unable to run in Discrimination Mode.\n";
	return (-1, -1);
    }
    elsif($n2 == 0){
	print LOG "\nNon-$category sequences have not been found in the dataset. Unable to run in Discrimination Mode.\n"; # Aborting...\n";
	print STDERR "\nNon-$category sequences have not been found in the dataset. Unable to run in Discrimination Mode.\n";
	return (-1, -1);
    }
    my $p1 = $n1/$n_tot;
    my $p2 = $n2/$n_tot;
    my $length = length($alignment[0]);
    my $i;
    my $j;
    my $aa;
    my $tam = scalar(@alignment);

    my @matrix = ();

    for($i = 0; $i < $tam; ++$i){
    	my @aux = split('', $alignment[$i]);
    	for($j = 0; $j < $length; ++$j){
            $matrix[$i][$j] = $aux[$j];
    	}
    }

    my $px;
    my $py;
    my $pi;
    my @values;
    my $aux;
    my %hash_alpa;
    my %hash_nalpa;
    my @values;
    my $mean = 0;
    my $dp = 0;

    for($i = 0; $i < $length; ++$i){
    	$pi = 0;
    	$px = 0;
    	$py = 0;
	%hash_alpa  = ();
	%hash_nalpa = ();
	my $aux_hash;
	if($file_type == 1){
	    $aux_hash = initialize_DNA();    
	}
	else{
    	    $aux_hash = initialize_protein();
	}
	%hash_alpa  = %{$aux_hash};
	%hash_nalpa = %{$aux_hash};
  	for($j = 0; $j < $tam; ++$j){
            if($j < $n1){
            	if(defined $hash_alpa{$matrix[$j][$i]}){
                    ++$hash_alpa{$matrix[$j][$i]};
            	}
            }
            else{
            	if(defined $hash_nalpa{$matrix[$j][$i]}){
                    ++$hash_nalpa{$matrix[$j][$i]};
            	}
            }
    	}
    	foreach my $aa (sort keys %hash_alpa){
            $aux = $hash_alpa{$aa};
            my $n = $hash_nalpa{$aa};
            my $tot = $aux + $n;
            if($aa eq '-'){
            	my $perc = ($aux/$n1)*100;
            	if($perc > $gap){
                    $pi = -1000;
                    last;
            	}
            }
	    elsif(lc($aa) eq 'x'){}
            $px = 0;
            $py = 0;
            if($tot > 0){
            	if($aux > 0){
                    my $frex = ($aux/$n_tot);
                    my $pa = ($tot/$n_tot);
                    my $palpa = $p1;
                    $px = ($frex) * log2($frex/($pa*$palpa));
            	}
            	if($n > 0){
                    my $frey = ($n/$n_tot);
                    my $pa = ($tot/$n_tot);
                    my $pnalpa = $p2;
                    $py = ($frey) * log2($frey/($pa*$pnalpa));
            	}
            }
            $pi += ($px+$py);
        }
        push (@values, $pi);
	if($pi != -1000){
	    $mean += $pi;
	}
    }
    $mean = $mean/$length;
    my $largest = largest_score(\@values);
    for($i = 0; $i < $length; ++$i){
	if($values[$i] != -1000){
	    $dp =+ ($values[$i] - $mean)**2;
	}
    }
    if(lc($normalize) eq "yes"){
	for($i = 0; $i < $length; ++$i){
	    if($values[$i] != -1000){
            	$values[$i] = $values[$i]/$largest;
	    }
    	}
    }
    $dp = $dp/($length-1);
    $dp = sqrt($dp);
    print LOG "\nYour category presents a mean MI score of $mean +- $dp\n";
    return \@values, \@gaps_only;	
}

####################################################################################
#         Calculates Sequence Harmony for each position of the alignment           #
####################################################################################

sub sh{
    my $aux1 = shift;
    my $aux2 = shift;
    my $gap  = shift;
    my $category = shift;

    my @aux_n = @{$aux1};
    my @aux_a = @{$aux2};
    my $t = scalar(@aux_n);
    my @cat1;
    my @seq1;
    my @cat2;
    my @seq2;
    my @entropy = ();
    for(my $i = 0; $i < scalar(@aux_n); ++$i){
        my $name = $aux_n[$i];
        $name =~ s/>//g;
        my @aux = split('_', $name);
        if(uc($aux[0]) eq uc($category)){
            push @cat1, $name;
            push @seq1, $aux_a[$i];
        }
        else{
            push @cat2, $name;
            push @seq2, $aux_a[$i];
        }
    }

    my @names = (@cat1, @cat2);

    my @gaps_only = @{get_gapOnly(\@seq1)};
    
    my @alignment = (@seq1, @seq2);

    my $n1 = scalar(@cat1);
    my $n2 = scalar(@cat2);
    if($n1 == 0){
	print LOG "$category sequences have not been found in the dataset. Unable to run in Discrimination Mode. Aborting...\n";
	close(LOG);
        die "$category sequences have not been found in the dataset. Unable to run in Discrimination Mode. Aborting...\n";
    }
    elsif($n2 == 0){
	print LOG "Non-$category sequences have not been found in the dataset. Unable to run in Discrimination Mode. Aborting...\n";
	close(LOG);
        die "Non-$category sequences have not been found in the dataset. Unable to run in Discrimination Mode. Aborting...\n";
    }
    my $length = length($alignment[0]);
    my $i;
    my $j;
    my $aa;
    my $tam = scalar(@alignment);

    my @matrix = ();

    for($i = 0; $i < $tam; ++$i){
        my @aux = split('', $alignment[$i]);
        for($j = 0; $j < $length; ++$j){
            $matrix[$i][$j] = $aux[$j];
        }
    }
   
    my %hash_a;
    my %hash_b;
    my $mean = 0;
    for($i = 0; $i < $length; ++$i){
        %hash_a  = ();
        %hash_b = ();
        my $aux_hash;
        if($file_type == 1){
            $aux_hash = initialize_DNA();
        }
        else{
            $aux_hash = initialize_protein();
        }
        %hash_a  = %{$aux_hash};
        %hash_b = %{$aux_hash};
        for($j = 0; $j < $tam; ++$j){
            if($j < $n1){
                if(defined $hash_a{$matrix[$j][$i]}){
                    ++$hash_a{$matrix[$j][$i]};
                }
            }
            else{
                if(defined $hash_b{$matrix[$j][$i]}){
                    ++$hash_b{$matrix[$j][$i]};
                }
            }
        }
	my $sa = 0;
        my $sb = 0;
        my $sab = 0;
        my $shi = 0;
        foreach my $key (sort keys %hash_a){
	    my $aux = $hash_a{$key};
            my $n = $hash_b{$key};
            my $tot = $aux + $n;
	    if(lc($key) eq 'x'){}
            elsif($key eq '-'){
                my $perc = ($aux/$n1)*100;
                if($perc > $gap){
		    $shi = -1000;                    
                    last;
                }
            }
            my $a = $hash_a{$key}/$n1;
            my $b = $hash_b{$key}/$n2;
            if($a != 0){
            	$sa += $a*log2($a); 
            }
            if($b != 0){
                $sb += $b*log2($b); 
            }
	    if($a+$b!=0){
                $sab += ($a+$b)*log2($a+$b);
            }
        }
	if($shi == -1000){
	    push @entropy, -1000;
	}
	else{
	    $sa *= -1;
            $sb *= -1;
            $sab *= -1;
            my $sh = (1/2)*($sab-$sa-$sb);
            $sh *= -1;
            $sh = 1 - $sh;
            $mean += $sh;
            push @entropy, $sh;
	}
    }

    my $dp = 0;
    $mean /= $length;
    for($i = 0; $i < scalar(@entropy); ++$i){
        $dp =+ ($entropy[$i] - $mean)**2;
    }
    $dp = sqrt($dp);
    print LOG "Your category presents a mean SH score of $mean +- $dp\n";
    
    return \@entropy, \@gaps_only;
}

######################################################################################
#                 Initialization of a protein frequency vector                       #
######################################################################################

sub initialize_protein{
    my %hash = (
           "A" => 0,
           "R" => 0,
           "N" => 0,
           "D" => 0,
           "C" => 0,
           "Q" => 0,
           "E" => 0,
           "G" => 0,
           "H" => 0,
           "I" => 0,
           "L" => 0,
           "K" => 0,
           "M" => 0,
           "F" => 0,
           "P" => 0,
           "S" => 0,
           "T" => 0,
           "W" => 0,
           "Y" => 0,
           "V" => 0,
           "-" => 0);

    return \%hash;
}

#####################################################################################
#                Initialization of a nucleotide frequency vector                    #
#####################################################################################

sub initialize_DNA{
    my %hash = (
           "A" => 0,
           "C" => 0,
           "G" => 0,
           "T" => 0,
           "-" => 0);

    return \%hash;
}

#####################################################################################
#                  Calculates the highest value in a set of values                  #
#####################################################################################

sub largest_score{
    my $aux = shift;
    my @array = @{$aux};
    my $largest = $array[0];
    for(my $i = 1; $i < scalar(@array); ++$i){
	if($array[$i] > $largest){
	    $largest = $array[$i];
	}
    }
    return $largest;
}

#####################################################################################
#                             Verifies if a value is unique                         #
#####################################################################################

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

#####################################################################################
#  Verifies the content of a input FASTA file (nucleotide or amino acid sequences)  #
#####################################################################################

sub verifyFileType{
    my $file = shift;
    my $count = 0;
    my $num_lines = 0;
    my $type = 1; # 1 - nucleotide; 2 - protein
    my $first = 0;
    open(my $vf, "$file");
    while(<$vf>){
        chomp($_);
        if($_ =~ /^>/){
	    if($first == 0){
		$first = 1;
	    }
	    else{
		if($count > 0){
                    $type = 2;
		}
		last;
	    }
	}
        else{
            if ((index($_, 'L') != -1) or (index($_, 'l') != -1)) {
                ++$count;
            }
            if ((index($_, 'V') != -1) or (index($_, 'v') != -1)) {
                ++$count;
            }
            if ((index($_, 'I') != -1) or (index($_, 'i') != -1)) {
                ++$count;
            }
            if ((index($_, 'P') != -1) or (index($_, 'p') != -1)) {
                ++$count;
            }
            if ((index($_, 'F') != -1) or (index($_, 'f') != -1)) {
                ++$count;
            }
            if ((index($_, 'S') != -1) or (index($_, 's') != -1)) {
                ++$count;
            }
            if ((index($_, 'Q') != -1) or (index($_, 'q') != -1)) {
                ++$count;
            }
            if ((index($_, 'D') != -1) or (index($_, 'd') != -1)) {
                ++$count;
            }
            if ((index($_, 'E') != -1) or (index($_, 'e') != -1)) {
                ++$count;
            }
            if ((index($_, 'R') != -1) or (index($_, 'r') != -1)) {
                ++$count;
            }
            if ((index($_, 'H') != -1) or (index($_, 'h') != -1)) {
                ++$count;
            }
            if ((index($_, 'M') != -1) or (index($_, 'm') != -1)) {
                ++$count;
            }
	    if ((index($_, 'X') != -1) or (index($_, 'x') != -1)) {
                ++$count;
            }
	    if ((index($_, 'N') != -1) or (index($_, 'n') != -1)) {
                ++$count;
            }
	    if ((index($_, 'K') != -1) or (index($_, 'k') != -1)) {
                ++$count;
            }
	    if ((index($_, 'P') != -1) or (index($_, 'p') != -1)) {
                ++$count;
            }
	    if ((index($_, 'T') != -1) or (index($_, 't') != -1)) {
                ++$count;
            }
	    if ((index($_, 'W') != -1) or (index($_, 'w') != -1)) {
                ++$count;
            }
	    if ((index($_, 'Y') != -1) or (index($_, 'y') != -1)) {
                ++$count;
            }
        }
    }
    close($file);
    if($count > 0){
	$type = 2;
    }
    return $type;
}

##########################################################################################
#        Prints the scores of each position of a alignment in a tabular file             #
#                                 (Conservation mode)                                    # 
##########################################################################################

sub printFileOnlyOption {
    my $out_file = shift;
    my $values = shift;
    my $selected = shift;
    open(FILE, ">$out_file");
    my @scores_v = @{$values};
    my @selected_v = @{$selected};
    print FILE "Position\tAll scores\tSelected scores\n";
    for($i = 0; $i < scalar(@scores_v); ++$i){
	my $p = $i+1;
    	if($scores_v[$i] == -1000){
            print FILE "$p\t0\t0\n";
    	}
        else{
            print FILE "$p\t$scores_v[$i]\t$selected_v[$i]\n";
    	}
    }	
    close(FILE);
}

##########################################################################################
#         Prints the scores of each position of a alignment in a tabular file            #
#                                (Discrimination mode)                                   #
##########################################################################################

sub printFileMOption {
    my $out_file = shift;
    my $values_mi_sh = shift;
    my $selected_mi = shift;
    my $selected_sh = shift;
    my $selected_mi_sh = shift;
    print LOG "Scores file: $out_file\n";
    open(FILE, ">$out_file");
    my @scores_v = @{$values_mi_sh};
    my @selected_mi_v = @{$selected_mi};
    my @selected_sh_v = @{$selected_sh};
    my @selected_mi_sh_v = @{$selected_mi_sh};
    print FILE "Position\tMI/SH average scores\tMI Selected scores\tSH selected scores\tMI/SH Selected scores\n";
    for($i = 0; $i < scalar(@scores_v); ++$i){
	my $pos = $i + 1;
  	if($scores_v[$i] == -1000){
            print FILE "$pos\t0\t$selected_mi_v[$i]\t$selected_sh_v[$i]\t$selected_mi_sh_v[$i]\n";
        }
        else{
            print FILE "$pos\t$scores_v[$i]\t$selected_mi_v[$i]\t$selected_sh_v[$i]\t$selected_mi_sh_v[$i]\n";
        }
    }
    close(FILE);
}

##########################################################################################
#           Selects valid blocks of the aligment based on the calculated scores          #
##########################################################################################

sub select_valid_blocks{
    my $aux_array = shift;
    my $thre = shift;
    my $win = shift;
    my $cat_name = shift;
    my $qnt = shift;
    my @v_scores = @{$aux_array};
    my @values_v = ();
    my $len_scores = scalar(@v_scores);
    for($i = 0; $i < $len_scores; ++$i){
    	push @values_v, $v_scores[$i];
    }

    my $len_values = scalar(@values_v);
    my $start;
    my $end;
    my $count;
    my $i, $j;
    my $gaps;
    my @valid_windows_m = ();
    for($i = 0; $i < $len_values; $i++){
    	$start = $i;
    	$end = ($start+$win) - 1;
    	$count = 0;
	$gaps = 0;
    	if($end <= $len_values - 1){
            for($j = $start; $j <= $end; $j++){
		if($values_v[$j] == -1000 and (lc($window_gaps eq "yes"))){
		    $gaps = 1;
		}	
            	if($values_v[$j] >= $thre){		
                    ++$count;
            	}
            }
	    if($gaps == 0){
		my $aux = ($count/$window_size)*100;
            	if(($count/$window_size)*100 >= $porc_threshold){
            	    $aux = $start."-".$end;
               	    push @valid_windows_m, $aux;
            	}
    	    }
    	}
	else{
	    last;
	}
    }

    $start = undef;
    $end   = undef;

    my @valid_blocks_m = ();
    foreach my $t (@valid_windows_m){
    	if((!defined $start) and (!defined $end)){
            @aux2 = split("-", $t);
            $start = $aux2[0];
            $end   = $aux2[1];
    	}
   	else{
            @aux2 = split("-", $t);
            $start_w = $aux2[0];
            $end_w   = $aux2[1];
            if(($start < $start_w and $start_w < $end) and ($end > $start_w and $end < $end_w)){
            	$end = $end_w;
            }
            else{
		my $length = ($end - $start) + 1;
            	if($length >= $min_size_block){
                    $aux = $start."-".$end;
                    push @valid_blocks_m, $aux;
            	}
            	$start = $start_w;
            	$end = $end_w;
           }
    	}
    }
    if(defined $start){
	my $length = ($end - $start) + 1;
        if($length >= $min_size_block){
    	    $aux = $start."-".$end;	
    	    push @valid_blocks_m, $aux;
	}
    }
    my $t = scalar(@valid_blocks_m);
    return (\@values_v, \@valid_blocks_m);
}

##########################################################################################
#           Selects valid blocks of the aligment based on the calculated scores          #
##########################################################################################

sub select_valid_blocks_disc{
    my $aux_array = shift;
    my $aux_gap_only = shift;
    my $thre = shift;
    my $win = shift;
    my $cat_name = shift;
    my $qnt = shift;
    my @v_scores = @{$aux_array};
    my @gaps_only = @{$aux_gap_only};
    my @values_v = ();
    my $len_scores = scalar(@v_scores);
    for($i = 0; $i < $len_scores; ++$i){
        push @values_v, $v_scores[$i];
    }

    # Select valid windows
    my $len_values = scalar(@values_v);
    my $start;
    my $end;
    my $count;
    my $i, $j;
    my $gaps;
    my @valid_windows_m = ();
    my $count_gaps;
    for($i = 0; $i < $len_values; $i++){
	if($gaps_only[$i] == 1){
	    next;
	}
        $start = $i;
        $end = ($start+$win) - 1;
        $count = 0;
        $gaps = 0;
        if($end <= $len_values - 1){
            my $count_size = 0;
            my $pos = $start;   	    
            while(($count_size < $win) and ($pos < ($len_values - 1) and $gaps == 0)){
                if(($values_v[$pos] == -1000) and (lc($window_gaps eq "yes"))){
                    $gaps = 1;
                }
                if($gaps_only[$pos] == 0){
                    ++$count_size;
                }
                if($gaps_only[$pos] == 0 and $values_v[$pos] >= $thre){
                    ++$count;
                }
                ++$pos;
            }
            $end = $pos-1;
            if($gaps == 0){
                my $aux = ($count/$win)*100;
		my $perc = ($count/$win)*100;
                if(($count/$win)*100 >= $porc_threshold){
                    $aux = $start."-".$end;
                    push @valid_windows_m, $aux;
                }
            }
        }
        else{
            last;
        }
    }
    $start = undef;
    $end   = undef;

    my @valid_blocks_m = ();
    foreach my $t (@valid_windows_m){
        if((!defined $start) and (!defined $end)){
            @aux2 = split("-", $t);
            $start = $aux2[0];
            $end   = $aux2[1];
        }
        else{
            @aux2 = split("-", $t);
            $start_w = $aux2[0];
            $end_w   = $aux2[1];
            if(($start < $start_w and $start_w < $end) and ($end > $start_w and $end < $end_w)){
                $end = $end_w;
            }
            else{
                my $length = ($end - $start) + 1;
                if($length >= $min_size_block){
                    $aux = $start."-".$end;
                    push @valid_blocks_m, $aux;
                }
                $start = $start_w;
                $end = $end_w;
           }
        }
    }
    if(defined $start){
        my $length = ($end - $start) + 1;
        if($length >= $min_size_block){
            $aux = $start."-".$end;
            push @valid_blocks_m, $aux;
        }
    }
    my $t = scalar(@valid_blocks_m);
    return (\@values_v, \@valid_blocks_m);
}

##########################################################################################
#              Joins alignment blocks found by the MI and SH methods                     #
##########################################################################################

sub blocks_union{
    my $mi = shift;
    my $sh = shift;
    my $int = shift;
    my $val = shift;
    my @values = @{$val};
    my @blocks_mi;
    my @blocks_sh;
    opendir(diretorio, "$mi");
    my @fastas_mi = readdir(diretorio);
    closedir(diretorio);

    opendir(diretorio, "$sh");
    my @fastas_sh = readdir(diretorio);
    closedir(diretorio);

    for(my $i = 0; $i < scalar(@fastas_mi); ++$i){
        if($fastas_mi[$i] =~ /(\d)_(\d+)-(\d+)/){
            my $str = $2."-".$3;
            push @blocks_mi, $str;
        }
    }
    for(my $i = 0; $i < scalar(@fastas_sh); ++$i){
        if($fastas_sh[$i] =~ /(\d)_(\d+)-(\d+)/){
            my $str = $2."-".$3;
            push @blocks_sh, $str;
        }
    }
    @blocks_mi = sort { $a <=> $b } @blocks_mi;
    @blocks_sh = sort { $a <=> $b } @blocks_sh;
    my $tm = scalar(@blocks_mi);
    my $ts = scalar(@blocks_sh);
    if($tm == 0){
	my @types = ();
	for($i = 0; $i < $ts; ++$i){
	    push @types, 2;
	}
	return (\@blocks_sh, \@types);	
    }
    elsif($ts == 0){
        my @types = ();
        for($i = 0; $i < $tm; ++$i){
            push @types, 1;
        }
        return (\@blocks_mi, \@types);
    }
    else{
    	my $i, $j;
    	my @aux = ();
    	my @remove = ();
    	for($i = 0; $i < scalar(@blocks_mi); ++$i){
            for($j = 0; $j < scalar(@blocks_sh); ++$j){
            	my @mi = split("-", $blocks_mi[$i]);
            	my $sm = $mi[0];
            	my $em = $mi[1];
            	my @sh = split("-", $blocks_sh[$j]);
            	my $ss = $sh[0];
            	my $es = $sh[1];
            	if($sm == $ss and $em == $es){
                    $blocks_mi[$i] = $blocks_mi[$i]."/1";
		    push @remove, $j;
                    last;
            	}
            	elsif($ss > $em){
                    $blocks_mi[$i] = $blocks_mi[$i]."/1";
                    last;
            	}
	        elsif(($sm <= $ss) and ($em >= $es) and ($es > $sm)){
		    $blocks_mi[$i] = $blocks_mi[$i]."/1";
                    push @remove, $j;
                    last;
	    	}
            }
    	}
    	my %params = map { $_ => 1 } @remove;
    	if(scalar(@blocks_sh) > 0){
            for($i = 0; $i < scalar(@blocks_sh); ++$i){
	    	if(defined $params{$i}){}
	    	else{
            	    $blocks_sh[$i] = $blocks_sh[$i]."/2";
	    	}
            }
    	}
    	my @union = ();
    	my @aux = (@blocks_mi,@blocks_sh);
    	@aux = sort { $a <=> $b } @aux;
    	my $length = scalar(@aux);
    	my $end  = undef;
    	my $last_s = undef;
    	my $last_e = undef;
    	my $last_t = undef;
    	for($i = 0; $i < $length; ++$i){	
            my @aux1 = split("-", $aux[$i]);
            my $s1 = $aux1[0];
            my @auxt1 = split("/", $aux1[1]);
            my $e1 = $auxt1[0];
            my $t1 = $auxt1[1];
	    if(!defined $last_s){
	    	$last_s = $s1;
	    	$last_e = $e1;
	    	$last_t = $t1;
	    	$end = $e1;
	    }
	    else{
	    	if((defined $end) and ($e1 <= $end)){           
            	}
            	elsif($s1 == $last_s and $e1 == $last_e){#equals
            	    next;
            	}
	    	elsif($last_s < $s1 and $last_e < $e1){#first before second
            	    my $len = $s1 - $last_e - 1;
            	    if((defined ($interval)) and ($len <= $int)){
		    	if(lc($window_gaps) eq "yes"){
		 	    my $gaps = 0;
                    	    for(my $j = ($last_e-1); $j < ($s1-1); ++$j){
                            	if($values[$j] > 0){
                            	    $gaps = 1;
                            	}
                    	    }
			    if($gaps == 0){
			    	$last_e = $e1;
                            	$last_t = 3;
                            	$end = $e1;
			    }
			    else{
			    	my $a = $last_s."-".$last_e."/".$last_t;
                   	    	push @union, $a;
                    	    	$last_s = $s1;
                    	    	$last_e = $e1;
                    	    	$last_t = $t1;
                    	    	$end = $e1;
			    }
		       	}
		    	else{
            	    	    $last_e = $e1;
            	   	    $last_t = 3;
            	    	    $end = $e1;
		    	}
            	    }
            	    else{
		    	my $a = $last_s."-".$last_e."/".$last_t;
                    	push @union, $a;
                    	$last_s = $s1;
            	    	$last_e = $e1;
            	    	$last_t = $t1;
            	    	$end = $e1;
                    	#next;
            	    }
            	}
	    	elsif($last_s <= $s1 and $s1 < $last_e and $last_e >= $e1){#Second inside first
            	}
	    	elsif($last_s <= $s1 and $s1 < $last_e and $last_e < $e1){#New coordinate
		    $last_s = $last_s;
                    $last_e = $e1;
                    $last_t = 3;
                    $end = $e1;
           	}
	    }
        }   
     	my $junction = $last_s."-".$last_e."/".$last_t;
     	push @union, $junction; 
     	@union = uniq(@union);
     	@union = sort { $a <=> $b } @union;
     	my @types = ();
     	for($i = 0; $i < scalar(@union); ++$i){
            my @aux = split("/", $union[$i]);
            push @types, $aux[1];
            $union[$i] = $aux[0];
            if($aux[1] eq '3'){
            }
       	}
     	return (\@union, \@types);
    }
}

##########################################################################################
#              Identifies block intersection between MI and SH methods                   #
##########################################################################################

sub setsIntersection{
    my $set1 = shift;
    my $set2 = shift;
    my @s1 = @{$set1};
    my @s2 = @{$set2};
    my @final = ();
    for(my $i = 0; $i < scalar(@s1); ++$i){
	my @aux = split("-", $s1[$i]);
	my $start1 = $aux[0];
	my $end1 = $aux[1];
	for(my $j = 0; $j < scalar(@s2); ++$j){
	    @aux = split("-", $s2[$j]);
            my $start2 = $aux[0];
            my $end2 = $aux[1];	    
	    if($start2 > $start1 and $start2 > $end1){
		next;
	    }
	    else{
		if($start1 == $start2 and $end1 == $end2){#equal
		    my $str = $start1."-".$end1;
		    push @final, $str;
   	    	}
    	        elsif($start1 <= $start2 and $start1 < $end2 and $end1 >= $end2){
		    my $str = $start2."-".$end2;
                    push @final, $str;    
    	    	}
    	    	elsif($start1 <= $start2 and $start2 < $end1 and $end1 < $end2){
		    my $str = $start1."-".$end1;
                    push @final, $str;
    	    	}
		elsif($start1 >= $start2 and $end1 >= $end2 and $start1 < $end2 ){
		    my $str = $start1."-".$end1;
                    push @final, $str;
		}
		elsif($start1 >= $start2 and $end1 < $end2 and $start1 < $end2 ){
                    my $str = $start1."-".$end1;
                    push @final, $str;
                }	
	    }
	}
    }
    return \@final;
}

##########################################################################################
#                        Removes redundant sequences in a block                          #
##########################################################################################

sub remove_redundancy{
    my $aux1 = shift;
    my $aux2 = shift;
    my $name = shift;
    my @nms = @{$aux1};
    my @alignments = @{$aux2};
    my @remove = (); # 0 - false; 1 - true
    my @ids = ();
    my $size = scalar(@nms);
    my $seq_length = length($alignment[0]);
    for(my $i = 0; $i < $size; ++$i){
    	$remove[$i] = 0
    }
    for(my $i = 0; $i < $size-1; ++$i){
    	if($remove[$i] == 1){
            next;
    	}
    	else{
            for(my $j = $i+1; $j < $size; ++$j){
            	my $result = calculate_identity($alignments[$i], $alignments[$j]);
            	if($result == 1){
		    $remove[$j] = 1;
                }
            }
    	}
    }
    my @aux_nms = ();
    my @alns = ();
    my $count = 0; 
    for(my $i = 0; $i < $size; ++$i){
	if($remove[$i] == 1){
	    if(defined $name){
		print LOG "Block $name: Sequence $nms[$i] discarded (redundant sequence)\n";
	    }
	    else{
	    	print LOG "Sequence $nms[$i] discarded (redundant sequence)\n";
	    }
	    ++$count;
	} 
 	else{
 	     push @aux_nms, $nms[$i];
	     push @alns, $alignments[$i];	
 	}
    }   
    return $count, \@aux_nms, \@alns;
}

##########################################################################################
#                          Removes redundant sequences in a MSA                          #
##########################################################################################

sub remove_redundancy_MSA{
    my $aux1 = shift;
    my $aux2 = shift;
    my @names = @{$aux1};
    my @alignment = @{$aux2};
    my $size = scalar(@names);
    my $seq_length = length($alignment[0]);
    my @aux_names = ();
    my @aux_seqs = ();
    my %hash = ();
    my @aux_names = ();
    my @aux_seqs = ();
    my @matrix = ();
    my $count_m = 0;
    my $count = 0;
    for(my $i = 0; $i < $size; ++$i){
     	if(!defined $hash{$alignment[$i]}){
            push @aux_names, $names[$i];
            push @aux_seqs, $alignment[$i];
            $hash{$alignment[$i]} = $names[$i];	
	    my @aux = split("", $alignment[$i]);
	    for(my $j = 0; $j < scalar(@aux); ++$j){
	    	$matrix[$count_m][$j] = $aux[$j];
	    }
	    ++$count_m;
     	}
	else{
	    ++$count;
	}
    }
    return $count, \@aux_names, \@aux_seqs, \@matrix;
}

##########################################################################################
#                Calculates identity percentage between two sequences                    #
##########################################################################################

sub calculate_identity{
    my $aux1 = shift;
    my $aux2 = shift;
    my @seq1 = split("", $aux1);
    my @seq2 = split("", $aux2);
    my $count = 0;
    my $len = scalar(@seq1);
    for(my $i = 0; $i < $len; $i++){
        if($seq1[$i] eq $seq2[$i]){
            ++$count;
        }
    }
    if($count == $len){
        return 1;
    }
    else{
        return 0;
    }
}

##########################################################################################
#                    Joins blocks within a maximum allowable distance                    #
##########################################################################################

sub interval_union{
    my $aux = shift;
    my $length = shift;
    my $aux_values = shift;
    my @blocks = @{$aux};
    my @values = @{$aux_values};
    my $start;
    my $end;
    my @aux_t;
    my @final_blocks = ();
    foreach my $t (@blocks){
    	if((!defined $start) and (!defined $end)){
            @aux2 = split("-", $t);
            $start = $aux2[0];
            $end   = $aux2[1];
            $count_t = 0;
            $tp = 0;
        }
        else{
            @aux2 = split("-", $t);
            $start_w = $aux2[0];
            $end_w   = $aux2[1];
            if(($start_w - $end - 1) <= $interval){
		if(lc($window_gaps) eq "yes"){
		    my $gaps = 0;
		    for(my $j = $end; $j < $start_w; ++$j){
			if($values[$j] == -1000){
			    $gaps = 1;
			}
		    }
		    if($gaps == 0){
			$end = $end_w;
		    }
		    else{
			$aux = $start."-".$end;
                	push @final_blocks, $aux;
                	$start = $start_w;
                	$end = $end_w;
		    }
		}
		else{ 
            	    $end = $end_w;
		}
            }
            else{
                $aux = $start."-".$end;
                push @final_blocks, $aux;
                $start = $start_w;
                $end = $end_w;
            }
        }
        ++$count_t;
    }
    $aux = $start."-".$end;
    push @final_blocks, $aux;
    return \@final_blocks;
}

##########################################################################################
#                            Removes ends of blocks with gaps                            #
##########################################################################################

sub trimming_zeros{
    my $start = shift;
    my $end = shift;
    my $aux_v = shift;
    my $min_length = shift;
    my @value = @{$aux_v};
    while($value[$start] == 0){
   	++$start;
    }
    while($value[$end] == 0){
        --$end;
    }
    return $start, $end;
}

##########################################################################################
#                                Removes gaps-only columns                               #
##########################################################################################

sub trimming_gaps_only{
    my $start = shift;
    my $end = shift;
    my $aux_v = shift;
    my $min_length = shift;
    my @gap = @{$aux_v};
    while($gap[$start] == 1){
        ++$start;
    }
    while($gap[$end] == 1){
        --$end;
    }
    return $start, $end;
}

##########################################################################################
#                                 Removes zero-score columns                             #
##########################################################################################

sub eliminate_zeros{
    my $aux_b = shift;
    my $aux_v = shift;
    my $min_length = shift;
    my @value = @{$aux_v};
    my @aux = @{$aux_b};
    my @new_blocks = ();
    foreach my $t (@aux){
	my @aux2 = split("-", $t);
        my $start = $aux2[0];
        my $end   = $aux2[1];
        while($value[$start] == 0 or $value[$start] == -1000){	   
            ++$start;
        }
        while($value[$end] == 0 or $value[$end] == -1000){
            --$end;
        }
	my $len = ($end - $start) + 1;
	if($len >= $min_length){
	    my $str = "$start-$end";
	    push @new_blocks, $str;
	}
    }
    return \@new_blocks;
}

##########################################################################################
#              Selects the best regions within a previously selected region              #
##########################################################################################

sub select_max_blocks{
    my $aux_b = shift;
    my $aux_v = shift;
    my $length = shift;
    my @final_blocks = ();
    my @value = @{$aux_v};
    print LOG "Selecting best subregions:\n";
    my @aux = @{$aux_b};
    foreach my $t (@aux){
        my @blocks;
        my @aux2 = split("-", $t);
        my $start = $aux2[0];
        my $end   = $aux2[1];
	my $as = $start + 1;
	my $ae = $end + 1;
	print LOG "Block: $as - $ae\n";
        my $size = $end - $start + 1;
        if($size > $length){
            for(my $i = $start; $i < $end; ++$i){
                my $s = $i;
                my $e = ($i+$length)-1;
                if($e > $end){
                    last;
                }
                my $sum = 0;
		for(my $j = $s; $j < $e; ++$j){
                    $sum += $value[$j];
                }
		my $string = "$sum-$s-$e";
		$as = $s + 1;
		$ae = $e + 1;
                print LOG "Region: $as - $ae: score sum=$sum\n";
                push @blocks, $string;
	    }
	    @aux2 = split("-", $blocks[0]);
            my $sum_max = $aux2[0];
            my $start_max = $aux2[1];
            my $end_max = $aux2[2];
            for(my $i = 1; $i < scalar(@blocks); ++$i){
                @aux2 = split("-", $blocks[$i]);
                if($aux2[0] > $sum_max){
                    $sum_max = $aux2[0];
                    $start_max = $aux2[1];
                    $end_max = $aux2[2];
                }
            }
	    $as = $start_max + 1;
	    $ae = $end_max + 1;
	    my ($a, $e) = trimming_zeros($start_max, $end_max, \@value, $min_size_block);
	    print LOG "Top score region: $as - $ae: score sum=$sum_max\n";	    
	    if(($e - $a + 1) >= $min_size_block){
            	my $tf = $a."-".$e;
            	push @final_blocks, $tf;
		++$e;
		++$a;
		print LOG "Selected region after zero-score trimming:$a-$e\n";
	    }
	    else{
		print LOG "Selected region is shorter than the minimum block size: discarded\n";
	    }
	    print LOG "\n";
        }
	else{
	    my $sum = 0;
	    for(my $j = $start; $j < $end; ++$j){
            	$sum += $value[$j];
            }
	    $as = $start + 1;
            $ae = $end + 1;
	    print LOG "Top score region: $as - $ae: score=$sum\n";
	    my ($a, $e) = trimming_zeros($start, $end, \@value, $min_size_block);
	    if(($e - $a + 1) >= $min_size_block){
                my $tf = $a."-".$e;
                push @final_blocks, $tf;
                ++$e;
                ++$a;
                print LOG "Selected region after zero-score trimming:$a-$e\n";
            }
            else{
                print LOG "Selected region is shorter than the minimum block size: discarded\n";
            }
	    print LOG "\n";
	}
    }
    return \@final_blocks;
}

##########################################################################################
#                       Creates FASTA files from selected blocks                         #
##########################################################################################

sub createOriginalBlocksFiles{
    my $dir = shift;
    my $aux_blocks = shift;
    my $aux_aln = shift;
    my $measure = shift;
    my $gap_perc = shift;
    my @blocks = @{$aux_blocks};
    my @alignments = @{$aux_aln};
    my $dir_fastas = $dir."/original_blocks";
    my $fasta_cat1;
    my $fasta_cat2;
    my $count = 1;
    my $len_aln = scalar(@alignments);
    foreach my $t (@blocks){
    	@aux2 = split("-", $t);
    	my $start = $aux2[0];
    	my $end   = $aux2[1];
    	my $as = $start + 1;
    	my $ae = $end + 1;
    	my $fasta = $dir."/".$outfile_name."_".$count."_".$as."-".$ae.".fasta";
	if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
            $fasta_cat1 = $dir."/".$outfile_name."_".$count."_".$as."-".$ae."_".$category.".fasta";
            $fasta_cat2 = $dir."/".$outfile_name."_".$count."_".$as."-".$ae."_non-".$category.".fasta";
            open(OUT2, ">$fasta_cat1");
            open(OUT3, ">$fasta_cat2");
   	}
	open(OUT, ">$fasta");
   	++$count;
	for($i = 0; $i < $len_aln; ++$i){
            my $len;
            my $seq;
            print OUT ">$names[$i]\n";
            $len = ($end-$start)+1;
            $seq = substr($alignment[$i], $start, $len);
            print OUT "$seq\n";
	    if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
            	my @aux = split('_', $names[$i]);
            	if(!(uc($aux[0]) eq uc($category))){
                    print OUT3 ">$names[$i]\n";
                    print OUT3 "$seq\n";
            	}
            	else{
                    print OUT2 ">$names[$i]\n";
                    print OUT2 "$seq\n";
            	}
             }
    	}
        close(OUT);
	if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
            close(OUT2);
            close(OUT3);
            if(-z $fasta_cat2){
            	system "rm $fasta_cat2";
            }
            if(-z $fasta_cat1){
            	system "rm $fasta_cat1";
            }
    	}
    	for($i = $start; $i <= $end; ++$i){
            $valid_blocks[$i] = $values[$i];
    	}
    }
    return \@valid_blocks;
}

##########################################################################################
#              Creates FASTA files from the selected blocks after cleaning               #
##########################################################################################

sub createFinalBlocksByOriginal{
    my $original_dir = shift;
    my $withoutRed_dir = shift;
    my $aux_values = shift;
    my $min_size_block = shift;
    my @values = @{$aux_values};
    opendir(diretorio, "$original_dir");
    my @fastas = readdir(diretorio);
    closedir(diretorio);
    my $size_selected = scalar(@fastas);
    if($size_selected == 0){
   	print STDERR "No block selected!\n";
   	print LOG "No block selected!\n";
   	my $file = $output."/scores.csv";
   	open(FILE, ">$file");
   	print FILE "Position\tAll scores\n";
    	for($i = 0; $i < scalar(@values); ++$i){
            my $p = $i+1;
            if($values[$i] == -1000){
            	print FILE "$p\t0\n";
            }
            else{
            	print FILE "$p\t$values[$i]\n";
            }
    	}
   	close(FILE);
   	close(LOG);
   	exit;
    }
    foreach my $file (@fastas){
    	if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
            if(!($file =~ /_$category.fasta/) or ($file =~ /non-/)){
            	next;
            }
    	}
	if(($file eq '.') or ($file eq '..') or (-z $file)){}
    	else{
  	     my $pos = rindex($file, ".");
             my $prefix = substr $file, 0, ($pos);
             # Aligning sequences to calculate their distance matrix
             $file = $original_dir."/".$file;
             my $start = 0;
             my $seq = "";
             my $name;
             my @seqs = ();
             my @nms = ();
	     open(FILE, "$file");
             while(<FILE>){
           	chomp($_);
            	if($_ =~ />/){
                    if($start == 0){
                    	$start = 1;
                    	$name = $_;
                    	$name =~ s/>//g;
                    	$name =~ s/^\s+|\s+$//g;
                     }
                     else{
                    	push @nms, $name;
                    	push @seqs, $seq;
                    	$seq = "";
                    	$name = $_;
                    	$name =~ s/>//g;
                    	$name =~ s/^\s+|\s+$//g;
                     }
            	}
            	else{
                    $seq .= $_;
            	}
            }
	    push @nms, $name;
            push @seqs, $seq;
            close(FILE);
	    my $t = scalar(@seqs);
	    my $prefix_name;
	    my $aux_s;
	    my $aux_e;
            my $num_seq;
	    if($prefix =~ /(\w+)_(\d+)_(\d+)-(\d+)/){
		$num_seq = $2;
            	$aux_s = $3;
                $aux_e = $4;
            }
	    $prefix_name = $aux_s."-".$aux_e;
	    my $file_name = $num_seq."_".$aux_s."-".$aux_e;
	    if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
		$file_name = $category."_".$file_name.".fasta";
	    }
	    else{
		if($measure eq "c"){
		    my $pos = index($input_file_prefix, "_");
                    if($pos == -1){
			$file_name =  $input_file_prefix."_".$file_name.".fasta";
                    }
                    else{
                        my $name = substr ($input_file_prefix, 0, $pos);
			$file_name =  $name."_".$file_name.".fasta";
                    }
                }		
	    }
            my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@nms, \@seqs, $prefix_name);
            @nms = @{$aux_name};
            @seqs = @{$aux_seqs};
	    my $seq_without_red = $withoutRed_dir."/".$file_name;
            open(OUT, ">$seq_without_red");
            my $i;
            for($i = 0; $i < scalar(@seqs); ++$i){
            	print OUT ">$nms[$i]\n$seqs[$i]\n";
            }
            close(OUT);
            # Aligning non-redundant sequences
            my $nseq = scalar(@nms); 
            my $min_nseq = 2;
            if(defined $amount){
            	$min_nseq = $amount;
            }	
            if($nseq < $min_nseq){
		print LOG "Number of sequences lower than $min_nseq: discarded \n\n";
		system "rm $seq_without_red";
            	next;
            }
	    else{
		my $aux = $withoutRed_dir."/".$prefix."_without_red_aln.fasta";
	        my ($aux1, $aux2) = read_fasta_alignment_short($seq_without_red);
           	my @names = @{$aux1};
           	my @alignment = @{$aux2};
           	my $size = scalar(@alignment);
           	my $len = length($alignment[0]);
    		my @del = ();
    		for(my $k = 0; $k < $len; ++$k){
        	    my @col = @{get_column($k, \@alignment)};
        	    my $gap = 100*(gap_percentage(\@col));
        	    if($gap == 100){
            		push @del, $k;
        	    }
    		}
    		for(my $k = 0; $k < $size; ++$k){
        	    my @aux = split("", $alignment[$k]);
        	    my @aux2 = ();
        	    for(my $d = 0; $d < scalar(@aux); ++$d){
            	    	if(contains($d, \@del) == 0){
                	    push @aux2, $aux[$d];
            		}
        	    }
        	    $alignment[$k] = join("", @aux2);
    		}
		$len = length($alignment[0]);		
           	if($len < $min_size_block){		    
		    system "rm $seq_without_red";
		    next;
           	}
           	#open(OUT, ">$aux");
		if(defined $gap_perc){
                    $len = length($alignment[0]);
		    my @aux_names = ();
		    my @aux_seqs = ();
		    for(my $j = 0; $j < $size; ++$j){
                	my @aux = split("", $alignment[$j]);
                	my $count = 0;
			my $count_x = 0;
                	for(my $k = 0; $k < $len; ++$k){
                    	    if($aux[$k] eq '-'){
                            	++$count;
                    	    }
			    elsif(lc($aux[$k]) eq 'x'){
                                ++$count_x;
                            }
                        }
                	my $perc = ($count/$len)*100;
			my $perc_x = ($count_x/$len)*100;
                	if($perc <= $gap_perc and $perc_x <= $gap_perc){
			    push @aux_names, $names[$j];
			    push @aux_seqs, $alignment[$j]; 
                    	    print OUT ">$names[$j]\n";
                    	    print OUT "$alignment[$j]\n";
                	}
                	else{
			    my $as;
			    my $ae;
			    if($prefix =~ /(\w+)_(\d+)-(\d+)/){
			    	$as = $2;
				$ae = $3;
			    }
			    if($perc > $gap_perc and $perc_x > $gap_perc){
                            	my $auxp = sprintf("%.2f", $perc);
                            	my $auxx = sprintf("%.2f", $perc_x);
				print LOG "Block $as - $ae: Sequence $names[$j] removed ($auxp\% of gaps in sequence and $auxx\% of X\'s in sequence)\n";
                            	next;
                            }
			    elsif($perc > $gap_perc){
                            	my $auxp = sprintf("%.2f", $perc);
                            	print LOG "Block $as - $ae: Sequence $names[$j] removed ($auxp\% of gaps in sequence)\n";
                            	next;
                            }					
			    elsif($perc_x > $gap_perc){
                            	my $auxx = sprintf("%.2f", $perc_x);
                            	print LOG "Block $as - $ae: Sequence $names[$j] removed  ($auxx\% of X\'s in sequence)\n";
                            	next;
                            }
                	}
		    }
		    my $nseq = scalar(@aux_names); 
                    my $min_nseq = 2;
                    if(defined $amount){
                        $min_nseq = $amount;
                    }
                    if($nseq < $min_nseq){
                        print LOG "Number of sequences lower than $min_nseq: discarded \n\n";
			system "rm $seq_without_red";
                        next;
                    }
                    else{
                        open(FILE, ">$aux");
                        for(my $k = 0; $k < scalar(@aux_names); ++$k){
                            print FILE ">$aux_names[$k]\n";
                            print FILE "$aux_seqs[$k]\n";
                        }
                        close(FILE);
			system "mv $aux $seq_without_red";
                    }
            	}
            	else{
		    open(OUT, ">$aux");
           	    for(my $j = 0; $j < $size; ++$j){
                    	print OUT ">$names[$j]\n$alignment[$j]\n";
            	    }
		    close(OUT);
                    system "mv $aux $seq_without_red";
		}
	    }
	}
    }
}


##########################################################################################
#                                Generates profile HMMs                                  #
##########################################################################################

sub createHMMs{
    my $aln_dir = shift;
    my $hmms_dir = shift;
    my $measure = shift;
    my $category = shift;
    my @final_coord = ();
    opendir(diretorio, "$aln_dir");
    my $msg = undef;
    my @fastas = readdir(diretorio);
    closedir(diretorio);
    foreach my $file (@fastas){
	if(($file eq '.') or ($file eq '..') or (-z $file)){}
	else{
	    my $pos = rindex($file, ".");
            my $prefix = substr $file, 0, ($pos);
	    my $caracter = substr $prefix, -1;
    	    if($caracter eq '_'){
        	my $l = length($prefix);
        	$prefix = substr $prefix, 0, ($l-1);
      	    }
	    $pos = index($prefix, "_");
	    my $hmm = $prefix.".hmm";
	    $hmm =~ s/\.fasta//gc;
	    $hmm = $hmms_dir."/".$hmm;
	    $file = $aln_dir."/".$file;
	    $hmmbuild =~ s/\"//g;
	    my $type;
	    if($file_type == 1){
		$type == "--dna";
	    }
	    else{
		$type = "--amino";
	    }
            my $resp = system "hmmbuild $type $hmmbuild $hmm $file > /dev/null";
	    if($resp > 0){
	        $msg = "ERROR: Could not run hmmbuild: hmmbuild $hmmbuild $hmm $file\n";
	    }
	    my $aux_name = $hmm;
            $aux_name =~ s/$hmms_dir//g;
            $aux_name =~ s/\///g;
            $aux_name =~ s/.hmm//g;
	    if($aux_name =~ /(\d+)_(\d+)-(\d+)/){
                my $str = $2."-".$3;
                push @final_coord, $str;
            }
	}
    } 
    return (\@final_coord, $msg); 
}

##########################################################################################
#                    Profile HMMs validation (Discrimination mode)                       #
##########################################################################################

sub hmmValidation{
    my $hmm_dir = shift;
    my $category = shift;
    my $num_cat = shift;
    my $seqs = shift;
    my $seq_type = shift;
    opendir(diretorio, "$hmm_dir");
    my @hmms = readdir(diretorio);
    closedir(diretorio);    
    my $analysis = $hmm_dir."/validation";
    system "mkdir $analysis";
    system "cp $seqs $analysis";
    $len = length($seqs);
    $seqs = substr $seqs, (rindex($seqs, "/") + 1), $len;
    $seqs = $analysis."/".$seqs;
    my $valid = $hmm_dir."/valid_HMMs";
    system "mkdir $valid";
    my $excluded = $hmm_dir."/excluded_HMMs";
    system "mkdir $excluded";
    print LOG "\nValidation of generated HMMs: \n\n";
    my $validation_file = $analysis."/results_".$category.".csv";
    open(VAL, ">$validation_file");
    print VAL "Model\t# $category\t% $category\n";
    foreach my $hmm (sort @hmms){
	if($hmm eq "." or $hmm eq ".."){
	    next;
	}
	print LOG "HMM: $hmm\n";
	my $pos = rindex($hmm, ".");
    	my $pos_f = rindex($hmm, "/");
	my $input_file_prefix; 
    	if($pos_f > -1){
            $pos -= ($pos_f+1);
            $input_file_prefix = substr $hmm, ($pos_f+1), ($pos);
    	}
    	else{
            $input_file_prefix = substr $hmm, 0, ($pos);
    	}
	my $analysis_dir = $analysis."/".$input_file_prefix;
	system "mkdir $analysis_dir";	
	my $file = $hmm_dir."/".$hmm;
	my $command;
        if($seq_type == 1){
	     $command = "nhmmer -T 1 --tblout $analysis_dir/results.tab -o $analysis_dir/results.txt $file $seqs 2>> /dev/null";
	}
	else{
	    $command = "hmmsearch -T 1 --tblout $analysis_dir/results.tab -o $analysis_dir/results.txt $file $seqs 2>> /dev/null";
	}
	system $command;	
	my $min_score = `grep "LENG" $file`;
	$min_score =~ s/LENG//;
	$min_score =~ s/\s+//g;
	$min_score *= $pontuation;
	open(FILE, "$analysis_dir/results.tab");
	my $cat_count = 0;
	my $cat_score;
	my $cat_evalue;
	my $other_score = undef;
	my $other_evalue = undef;
	my $intersec;
	my $first = 0;
	my $discart = 0;
	while(<FILE>){
	    chomp($_);
	    if($_ =~ /#/){}	    
	    else{
		my @aux = split(" ", $_);
		my $name = $aux[0];
		$name =~ s/>//g;
            	my @aux_name = split('_', $name);
		if($first == 0 and !(lc($aux_name[0]) eq lc($category))){
		    $discart = 1;
		    last;
		}
		++$first;
		if(lc($aux_name[0]) eq lc($category)){
		    ++$cat_count;
		    if($seq_type == 1){
			$cat_score = $aux[13];
                    	$cat_evalue = $aux[12];
		    }
		    else{
		    	$cat_score = $aux[8];
		    	$cat_evalue = $aux[7];
		    }
	    	}
	        else{		
		    if($seq_type == 1){
                    	$other_score = $aux[13];
                        $other_evalue = $aux[12];
                    }
                    else{
		       	$other_score = $aux[8];
                    	$other_evalue = $aux[7];
		    }
		    last;
	    	}
	    }
	}
	close(FILE);
	my $aux_thre = $threshold_validation/100; 
	my $threshold = $cat_score*$aux_thre;
        my $evalue_thre = $cat_evalue*$aux_thre;
        print LOG "\tNumber of sequences of group $category in the MSA: $num_cat\n";
	my $aux_file;
	if($discart == 1){
	     $aux_file = $excluded."/".$hmm;
             print LOG "\tStatus: Discarded - Top score hit does not belong to the selected group ($category)\n\n";
	     print VAL "$input_file_prefix\t0\t0.0\n";
	}
	elsif(defined $other_score){
	    my $delta = ($cat_score - $other_score)*0.8;
	    $threshold = $other_score + $delta;
	    if(($min_score > $threshold) and ($full_length eq "no")){
	    	print LOG "\tSuggested score is lower ($threshold) than minimum score ($min_score). New suggested score = $min_score\n";
		$threshold = $min_score;
	    }
	    my $aux = ($other_score*100/$cat_score);
	    if($aux <= $threshold_validation){
	        my $count = 0;
		open(FILE, "$analysis_dir/results.tab");
		while(<FILE>){
		    chomp($_);
 	            if($_ =~ /#/){}
            	    else{
                    	my @aux = split(" ", $_);
			my $name = $aux[0];
	                $name =~ s/>//g;
        	        my @aux_name = split('_', $name);
			if(lc($aux_name[0]) eq lc($category)){
                    	    if($seq_type == 1){
                            	if($aux[13] >= $threshold){
				    ++$count;
                    		}
			    }
                    	    else{
                            	if($aux[8] >= $threshold){
			   	    ++$count;
                    		}
		    	    }
		   	}
		    }
		}
		my $perc = ($count/$num_cat)*100;
   	        my $print_perc = sprintf("%.2f",$perc);
        	print VAL "$input_file_prefix\t$count\t$print_perc\n";
		print LOG "\tNumber of sequences of group $category detected by the model: $count\n";
        	print LOG "\tLowest score of group $category: $cat_score\n";
		print LOG "\tHighest score of the remaining sequences: $other_score\n";
                print LOG "\tRecommended cutoff score: $threshold\n";
		close(FILE);
		if(($count*100/$num_cat) >= $amount_detection){
		    $aux_file = $valid."/".$hmm;
                    print LOG "\tStatus: Valid\n\n";
		}
		else{
	            print LOG "\tStatus: Discarded - model did not detect >=$amount_detection% of the category $category of the training set\n\n";
        	    $aux_file = $excluded."/".$hmm;
		}
	    }
	    else{
	    	$aux_file = $excluded."/".$hmm;
		open(FILE, "$analysis_dir/results.tab");
		$count = 0;
                while(<FILE>){
                    chomp($_);
                    if($_ =~ /#/){}
                    else{
                    	my @aux = split(" ", $_);
                        my $name = $aux[0];
                        $name =~ s/>//g;
                        my @aux_name = split('_', $name);
                        if(lc($aux_name[0]) eq lc($category)){
                            if($seq_type == 1){
                            	if($aux[13] >= $threshold){
                                    ++$count;
                                }
                            }
                            else{
                            	if($aux[8] >= $threshold){
                                     ++$count;
                                }
                            }
                        }
                   }
               }
	        print LOG "\tNumber of sequences of group $category detected by the model: $count\n";
                print LOG "\tLowest score of group $category: $cat_score\n";
                print LOG "\tHighest score of the remaining sequences: $other_score\n";
                print LOG "\tRecommended cutoff score: $threshold\n";
                print LOG "\tStatus: Discarded - the highest score observed of non-selected group is >=$threshold_validation% of the lowest score of group $category\n\n";
                my $perc = ($count/$num_cat)*100;
                my $print_perc = sprintf("%.2f",$perc);
                print VAL "$input_file_prefix\t$count\t$print_perc\n";
	    }		
	}
	else{
	    if($min_score > $threshold){
	        print LOG "\tSuggested score is lower ($threshold) than minimum score ($min_score). New suggested score = $min_score\n";
                $threshold = $min_score;
            }
	    my $count = 0;
            open(FILE, "$analysis_dir/results.tab");
	    while(<FILE>){
            	chomp($_);
                if($_ =~ /#/){}
                else{
                    my @aux = split(" ", $_);
                    my $name = $aux[0];
                    $name =~ s/>//g;
                    my @aux_name = split('_', $name);
                    if(lc($aux_name[0]) eq lc($category)){
                    	if($seq_type == 1){
                            if($aux[13] >= $threshold){
                               ++$count;
                            }
                        }
                        else{
                            if($aux[8] >= $threshold){
                                ++$count;
                            }
                        }
                    }
                }
	    }
	    my $perc = ($count/$num_cat)*100;
            my $print_perc = sprintf("%.2f",$perc);
            print VAL "$input_file_prefix\t$count\t$print_perc\n";
	    print LOG "\tNumber of sequences of group $category detected by the model: $count\n";
            print LOG "\tLowest score of group $category: $cat_score\n";
            print LOG "\tRecommended cutoff score: $threshold\n";
            close(FILE);
            if(($count*100/$num_cat) >= $amount_detection){
            	$aux_file = $valid."/".$hmm;
                print LOG "\tStatus: Valid\n\n";
            }
	    else{
	    	print LOG "\tStatus: Discarded - model did not detect >=$amount_detection% of the category $category of the training set\n\n";
                $aux_file = $excluded."/".$hmm;
	    }
	}
  	if($cutoff_score eq "yes"){	
	    open(FILE, "$file");
	    open(OUT, ">$aux_file");
	    while(<FILE>){
     	    	chomp($_);
     	    	if($_ =~ /STATS LOCAL FORWARD/){
		    my $string = $_."\nCUTOFF SCORE\t$threshold";
            	    print OUT "$string\n";
      	    	}
            	else{
            	    print OUT "$_\n";
      	    	}
	    }
	    close(FILE);
	    close(OUT);
	    system "rm $file";
	}
	else{
	    system "mv $file $aux_file";
	}
    }
    opendir(diretorio, "$valid");
    @hmms = readdir(diretorio);
    closedir(diretorio);
    my @regions = ();
    foreach my $hmm (@hmms){
	 my $aux = $outfile_name."_";
	 my $aux_name = $hmm;
         $aux_name =~ s/$aux//g;
         if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
             $aux = "_".$category;
             $aux_name =~ s/$aux//g;
         }
         if($aux_name =~ /(\d)_(\d+)-(\d+)/){
             my $str = $2."-".$3;
             push @regions, $str;
         }
    }
    close(VAL);
    if($clean eq "yes"){
	system "rm -rf $excluded";
    }
    return \@regions;
}

##########################################################################################
#                      Profile HMMs validation (Conservation mode)                       #
##########################################################################################

sub hmmValidation_conservation{
    my $hmm_dir = shift;
    my $num_cat = shift;
    my $seqs = shift;
    my $seq_type = shift;
    my $categ = shift;
    my $am = shift;
    my @cats; 
    my @amount;
    if(defined $categ){
	@cats = @{$categ};
	@amount = @{$am};
    }
    opendir(diretorio, "$hmm_dir");
    my @hmms = readdir(diretorio);
    closedir(diretorio);
    my $analysis = $hmm_dir."/validation";
    system "mkdir $analysis";
    system "cp $seqs $analysis";
    $len = length($seqs);
    $seqs = substr $seqs, (rindex($seqs, "/") + 1), $len;
    $seqs = $analysis."/".$seqs;
    my $validation_file = $analysis."/results.csv";
    open(VAL, ">$validation_file");
    print VAL "Model";
    if(defined $categ){
    	for(my $i = 0; $i < scalar(@cats); ++$i){
             print VAL "\t#$cats[$i]\t%$cats[$i]";
    	}
    }
    print VAL "\tTotal # of sequences\t% of all sequences\n";
    my $valid = $hmm_dir."/valid_HMMs";
    system "mkdir $valid";
    my $excluded = $hmm_dir."/excluded_HMMs";
    system "mkdir $excluded";
    print LOG "\nValidation of generated HMMs: \n\n";
    foreach my $hmm (sort @hmms){
        if($hmm eq "." or $hmm eq ".."){
            next;
        }
        print LOG "HMM: $hmm\n";
        my $pos = rindex($hmm, ".");
        my $pos_f = rindex($hmm, "/");
        my $input_file_prefix;
        if($pos_f > -1){
            $pos -= ($pos_f+1);
            $input_file_prefix = substr $hmm, ($pos_f+1), ($pos);
        }
        else{
            $input_file_prefix = substr $hmm, 0, ($pos);
        }
        my $cat_score;
        my $cat_evalue;
        my $analysis_dir = $analysis."/".$input_file_prefix;
        system "mkdir $analysis_dir";
	my $validation_file;
        my $file = $hmm_dir."/".$hmm;
   	my $command;
        my $threshold;
        if($seq_type == 1){
	    $command = "nhmmer -T 1 --tblout $analysis_dir/results.tab -o $analysis_dir/results.txt $file $seqs ";
	}
	else{
            $command = "hmmsearch -T 1 --tblout $analysis_dir/results.tab -o $analysis_dir/results.txt $file $seqs 2>> /dev/null";
	}
        system $command;
	open(FILE, "$analysis_dir/results.tab");
	my $aux_file;
	my $valid_hmm = 0;
     	my @names = ();
	my @scores = ();
	my $cat_count = 0;	
        my $min_score = `grep "LENG" $file`;
        $min_score =~ s/LENG//;
        $min_score =~ s/\s+//g;
        $min_score *= $pontuation;
  	while(<FILE>){
            chomp($_);
            if($_ =~ /#/){}
            else{
                my @aux = split(" ", $_);
		push @names, $aux[0];
		if($seq_type == 1){
                    $cat_score = $aux[13];
                    $cat_evalue = $aux[12];
		    push @scores, $aux[13];
                }
	        else{
		    $cat_score = $aux[8];
                    $cat_evalue = $aux[7];
		    push @scores, $aux[8];
		}
		++$cat_count;
	    }
	}
	close(FILE);
	$threshold = $cat_score * 0.8;
	if($threshold < $min_score){# and ($full_length eq "no")){
	    print LOG "\tSuggested score is lower ($threshold) than minimum score ($min_score). New suggested score = $min_score\n";
	    $threshold = $min_score;	    
        }
	my $cats_name = "";
	if(defined $categ){
	    print VAL "$input_file_prefix\t";
	    my %aux_hash = ();
	    for(my $i = 0; $i < scalar(@cats); ++$i){		
	    	my $count = 0;
	    	my $name = lc($cats[$i]);
	    	for(my $j = 0; $j < scalar(@names); ++$j){
		    if((lc($names[$j]) =~ /^$name/) and ($scores[$j] >= $threshold)){
		    	++$count;
		    }	
	    	}
		my $perc;
		if($amount[$i] == 0){
		    $perc = 0;
		}
		else{
		    $perc = ($count/$amount[$i])*100;
		}
		my $formated = sprintf("%.2f",$perc);
		print VAL "$count\t$formated\t";
	    	if($perc < $amount_detection){
		    $valid_hmm = 1;
		    $cats_name .= $cats[$i]."\t";
	    	}
	    }
        }	    
        print LOG "\tNumber of sequences in the MSA: $num_cat\n";
	my $continue = 0;
	    if(defined $cat){
		if($valid_hmm == 1){
		    print LOG "\tRecommended cutoff score: $threshold\n";
		    my @aux = split("\t", $cats_name);
		    my $tot_cat = scalar(@cats);
		    my $not_detect = scalar(@aux);
		    my $perc_cat = (($tot_cat - $not_detect)/$tot_cat)*100;
		    if($perc_cat < $category_percentage){
			print LOG "\tStatus: Discarded - model did not detect >= $category_percentage of categories\n";
		    	if(scalar(@aux) == 1){
			    print LOG "\t	The model did not detect >=$amount_detection% of category $aux[0] of the training set\n\n";
		    	}
		        else{
			    my $str = "";
			    my @aux = split("\t", $cats_name);
			    my $size = scalar(@aux);
			    for(my $j = 0; $j < $size-2; ++$j){
			    	$str .= $aux[$j].",";
			    }
			    $str = $aux[$size-2]." and  ".$aux[$size-1];
		    	    print LOG "\t	The model did not detect >=$amount_detection% of categories $str of the training set\n\n";
		    	}
			$continue = 1;
		    }
		    else{
			$continue = 0
		    }
                    $aux_file = $excluded."/".$hmm;
		    my $count = 0;
                    open(FILE, "$analysis_dir/results.tab");
                    while(<FILE>){
                     	chomp($_);
                     	if($_ =~ /#/){}
                     	else{
                            my @aux = split(" ", $_);
                            if($seq_type == 1){
                            	if($aux[13] >= $threshold){
                                    ++$count;
                            	}
                            }
                            else{
                            	if($aux[8] >= $threshold){
                                    ++$count;
                            	}
                            }
                    	}
                    }
		    close(FILE);
                    my $total_perc = ($count/$num_cat)*100;
                    my $formated = sprintf("%.1f",$total_perc);
                    print VAL "$count\t$formated\n";
		}
		else{
		    $continue = 0;
		}
	    }	    
	    if($continue == 0){
		my $count = 0;
                open(FILE, "$analysis_dir/results.tab");
                while(<FILE>){
                     chomp($_);
                     if($_ =~ /#/){}
                     else{
                     	my @aux = split(" ", $_);
                        if($seq_type == 1){
                            if($aux[13] >= $threshold){
                            	++$count;
                            }
                        }
                        else{
                            if($aux[8] >= $threshold){
                            	++$count;
                            }
                        }
                    }
                }
		my $total_perc = ($count/$num_cat)*100;
        	my $formated = sprintf("%.1f",$total_perc);
        	print VAL "$count\t$formated\n";
		print LOG "\tNumber of sequences detected by the model: $count\n";
                print LOG "\tRecommended cutoff score: $threshold\n";
                close(FILE);
		my $total_perc = ($count/$num_cat)*100;
		if( $total_perc >= $amount_detection){
		    print LOG "\tStatus: valid\n\n";
		    $aux_file = $valid."/".$hmm;
	    	}
	    	else{
		    print LOG "\tStatus: Discarded - model did not detect >=$amount_detection% of the total number of sequences of the training set\n\n";
                    $aux_file = $excluded."/".$hmm;
	    	}
	    }
        if($cutoff_score eq "yes"){
	    open(FILE, "$file");
            open(OUT, ">$aux_file");
            while(<FILE>){
            	chomp($_);
            	if($_ =~ /STATS LOCAL FORWARD/){
                    my $string = $_."\nCUTOFF SCORE\t$threshold";
                    print OUT "$string\n";
            	}
            	else{
                    print OUT "$_\n";
            	}
            }
            close(FILE);
            close(OUT);
            system "rm $file";
	}
	else{
	    system "mv $file $aux_file";
	}
    }
    if(defined $cat){
	close(VAL);
    }
    opendir(diretorio, "$valid");
    @hmms = readdir(diretorio);
    closedir(diretorio);
    my @regions = ();
    foreach my $hmm (@hmms){
         my $aux = $outfile_name."_";
         my $aux_name = $hmm;
         $aux_name =~ s/$aux//g;
         if($aux_name =~ /(\d)_(\d+)-(\d+)/){
             my $str = $2."-".$3;
             push @regions, $str;
         }
    }
    return \@regions;
}

##########################################################################################
#                Creates FASTA files from the best selected sub-regions                  #
##########################################################################################

sub createFinalBlocksSelectingBestRegion{
    my $original_dir = shift;
    my $aux_values = shift;
    my $min_size_block = shift;
    my $max_block_size = shift;
    my $measure = shift;
    my $category = shift;
    my $prefix = shift;
    my $withoutRed_dir = shift;
    my $category = shift;
    my $aux_gaps_only = shift;
    my @value = @{$aux_values};
    my @gaps_only = @{$aux_gaps_only}; 
    opendir(diretorio, "$original_dir");
    my @fastas = readdir(diretorio);
    closedir(diretorio);
    my $size_selected = scalar(@fastas);
    if($size_selected == 0){
        print STDERR "No block selected!\n";
        print LOG "No block selected!\n";
        my $file = $output."/scores.csv";
        open(FILE, ">$file");
        print FILE "Position\tAll scores\n";
        for($i = 0; $i < scalar(@values); ++$i){
            my $p = $i+1;
            if($values[$i] == -1000){
                print FILE "$p\t0\n";
            }
            else{
                print FILE "$p\t$values[$i]\n";
            }
        }
        close(FILE);
        close(LOG);
        exit;
    }

    @fastas = sort {$a <=> $b or $a cmp $b} @fastas; #sort { $a <=> $b } @fastas;
    my %hash = ();
    foreach my $file (@fastas){
        if((lc($measure) eq "mi") or (lc($measure) eq "sh") or (lc($measure) eq "d")){
            if(!($file =~ /_$category.fasta/) or ($file =~ /non-/)){
                next;
            }
        }
        if(($file eq '.') or ($file eq '..') or (-z $file)){}
        else{	   
	    my $aux_name = $file;
	    $aux_name =~ s/$prefix//gci;
	    $aux_name =~ s/^_//g; 
	    if(defined $category){
	    	$aux_name =~ s/_$category//gi;
	    }
	    $aux_name =~ s/\.fasta//gci;
	    $aux_name =~ s/$prefix//gci;
	    if($aux_name =~ /(\d+)_(\d+)-(\d+)/){
		$hash{$file} = $2."-".$3;
	    }	    
	}
    }
    my @final_coord = ();
    my $num_seq = 0;
    foreach my $key (sort {$hash{$a} <=> $hash{$b} or $a cmp $b}  keys %hash){
	++$num_seq;
	my @aux = split("-", $hash{$key});
	my $start = $aux[0] - 1;
	my $end = $aux[1] - 1;
	my $as = $start + 1;
	my $ae = $end + 1;
	print LOG "\nBlock: $as - $ae\n";
	my $file = $original_dir."/".$key;
	my ($aux1, $aux2) = read_fasta_alignment_short($file);
        my @names = @{$aux1};
        my @alignment = @{$aux2};
	my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@names, \@alignment, undef);
        @names = @{$aux_name};
        @alignment = @{$aux_seqs};
	my $nseq = scalar(@names); 
        my $min_nseq = 2;
        if(defined $amount){
            $min_nseq = $amount;
        }
        if($nseq < $min_nseq){
            print LOG "Number of sequences lower than $min_nseq: discarded\n";
            next;
        }
	
        my $size = scalar(@alignment);
        my $len = length($alignment[0]);
 	my @delete = ();
        for(my $i = 0; $i < $len; ++$i){
            my @col = @{get_column($i, \@alignment)};
            my $gap = 100*(gap_percentage(\@col));
            if($gap == 100){
		push @delete, ($start+$i);
            }
            --$i;
            --$len;
        }
	$len = length($alignment[0]);
        if($len < $min_size_block){
            next;
        }
 	my $size = $end - $start + 1;
	my $count = 0;
	my $len_del = scalar(@delete);
        if($len_del > 0){
            for(my $k = 0; $k > $len_del; ++$k){
            	if($delete[$k] >= $start and $delete[$k] <= $end){
                    $count++;
                }
            }
        }
        if($size > $max_block_size){
	    my ($as, $ae) = selectBestRegion($start, $end, \@value, $count, \@gaps_only);
            my ($a, $e) = trimming_zeros(($as-1), ($ae-1), \@value, $min_size_block);
	    # Removes redundant sequences
	    my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@names, \@alignment, undef);
            @names = @{$aux_name};
            @alignment = @{$aux_seqs};
	    my $nseq = scalar(@names); 
            my $min_nseq = 2;
            if(defined $amount){
            	$min_nseq = $amount;
            }
            if($nseq < $min_nseq){
            	print LOG "Number of sequences lower than $min_nseq: discarded\n";
            	next;
            }
	
	    if(($e - $a + 1) >= $min_size_block){
		my $aux_start = $a - $start;
		my $aux_end = $e - $start;
		++$a;
		++$e;
                print LOG "Selected region after zero-score trimming:$a-$e\n";
		my $file_new;
		if($measure eq "c"){
		    my $pos = index($input_file_prefix, "_");
		    if($pos == -1){
			$file_new =  $withoutRed_dir."/".$input_file_prefix."_".$num_seq."_".$a."-".$e.".fasta";
		    }
		    else{
			my $name = substr ($input_file_prefix, 0, $pos);
		        $file_new =  $withoutRed_dir."/".$name."_".$num_seq."_".$a."-".$e.".fasta";
		    }
		}
		else{
		    $file_new =  $withoutRed_dir."/".$category."_".$num_seq."_".$a."-".$e.".fasta";
		}
		my @new_seqs = @{selectSequences(\@alignment, $start, $end, $a, $e)};
           	$len = length($new_seqs[0]);
           	if($len < $min_size_block){
		    print LOG "Selected region is shorter than the minimum block size after removing gap-only columns: discarded\n";
                    next;
           	}
		my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@names, \@new_seqs, undef);
                my @new_names = @{$aux_name};
                @new_seqs = @{$aux_seqs};
		my $nseq = scalar(@names); 
        	my $min_nseq = 2;
        	if(defined $amount){
            	    $min_nseq = $amount;
        	}
        	if($nseq < $min_nseq){
            	    print LOG "Number of sequences lower than $min_nseq: discarded\n";
            	    next;
        	}
		if(defined $gap_perc){
		    my ($aux_names, $aux_alns) = calculateGapPercentage(\@new_names, \@new_seqs, $start, $end);
		    my @aux_n = @{$aux_names};
                    my @aux_s = @{$aux_alns};
		    my $min_nseq = 2;
                    if(defined $amount){
                        $min_nseq = $amount;
                    }
                    if(scalar(@aux_n) >= $min_nseq){
                        open(FILE, ">$file_new");
                        for(my $k = 0; $k < scalar(@aux_n); ++$k){
                            print FILE ">$aux_n[$k]\n";
                            print FILE "$aux_s[$k]\n";
                        }
                        close(FILE);
                    }
		}
            	else{
		    open(FILE, ">$file_new");
		    for(my $k = 0; $k < scalar(@new_names); ++$k){
		    	print FILE ">$new_names[$k]\n";
                    	print FILE "$new_seqs[$k]\n"; 
		    }
		    close(FILE);
                }
		if(lc($more_blocks) eq "yes"){
		    print LOG "Trying to select more subregions:\n";		
		    my @aux_blocks = @{searchForMoreRegions($start, $end, ($a - 1), ($e - 1), \@value, \@gaps_only)};		
		    if(scalar(@aux_blocks) == 0){
		    	print LOG "No new subregion found!\n\n";
		    }
		    else{
			my $last = $e-1;
			my $first = $a-1;
			@aux_blocks = @{sortCoordinates(\@aux_blocks)};
		    	for(my $l = 0; $l < scalar(@aux_blocks); ++$l){
		    	    my @aux_b = split("-", $aux_blocks[$l]);
		    	    my $start_b = $aux_b[0];
			    if(($start_b == $last) or ($start_b == $e - 2)){
				++$start_b;
			    }
			    elsif(($start_b < $last) and (($last - $start_b) < 5)){
				my $dif = $last-$start_b;
				$start_b += ($dif + 1);
			    }
		    	    my $end_b = $aux_b[1];
			    if(($end_b == ($a - 1)) or ($end_b == $first)){
				--$end_b;
			    }
			    if(($end_b - $start_b + 1) < $min_size_block){
				next;
			    }
			    my ($new_s, $new_e) = trimming_zeros($start_b, $end_b, \@value, $min_size_block);
			    if(($new_e - $new_s + 1) >= $min_size_block){
		    	    	my $aux_s = $new_s;
		    	    	my $aux_e = $new_e;
				my ($ts, $te) = trimming_gaps_only($aux_s, $aux_e, \@gaps_only, $min_size_block);
		    	    	print LOG "New block: $ts - $te\n";
		    	    	++$num_seq;
		    	    	my $resp = generateRegion($start, $end, $ts, $te, \@value, \@names, \@alignment, $withoutRed_dir, $num_seq, $category);
		    	    	if($resp == -1){
			    	    --$num_seq;
		    	    	}
			    	$last = $new_e;
			    	$first = $new_s;
			    }
			    else{
				print LOG "Selected region is shorter than the minimum block size: discarded\n";
			    }
			}
		    }
		}
            }
            else{
                print LOG "Selected region is shorter than the minimum block size: discarded\n";
            }
	}	
	else{
	    my $sum = 0;
            for(my $j = $start; $j < $end; ++$j){	
		if($value[$j] != -1000){
                    $sum += $value[$j];
		}
            }
            $as = $start + 1;
            $ae = $end + 1;
            print LOG "Top score region: $as - $ae: score=$sum\n";
            my ($a, $e) = trimming_zeros($start, $end, \@value, $min_size_block);
            if(($e - $a + 1) >= $min_size_block){
		my $aux_start = 0;
                my $aux_end = 0;
                if($a > $start){ 
                    $aux_start = $a - $start;
                }
                if($e < $end){ 
                    $aux_end = ($end - $e);
                }
		$a++;
		$e++;
                print LOG "Selected region after zero-score trimming:$a-$e\n";
		my $file_new;
		if($measure eq "c"){
                    if($input_file_prefix =~ /_/){
                        my @aux = split("_", $input_file_prefix);
                        $file_new =  $withoutRed_dir."/".$aux[0]."_".$num_seq."_".$a."-".$e.".fasta";
                    }
                    else{
                        $file_new =  $withoutRed_dir."/".$category."_".$num_seq."_".$a."-".$e.".fasta";
                    }
                }
                else{
			$file_new =  $withoutRed_dir."/".$category."_".$num_seq."_".$a."-".$e.".fasta";
		}
		my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@names, \@alignment, undef);
		@names = @{$aux_name};
		@alignment = @{$aux_seqs};
		my $nseq = scalar(@names); 
        	my $min_nseq = 2;
        	if(defined $amount){
            	    $min_nseq = $amount;
        	}
        	if($nseq < $min_nseq){
            	     print LOG "Number of sequences lower than $min_nseq: discarded\n";
            	     next;
        	}

		my @new_seqs = ();
                for(my $k = 0; $k < scalar(@names); ++$k){
		    my $size = length($alignment[$k]);
                    my $len = $size - $aux_end - $aux_start + 1;		    
                    my $str = substr $alignment[$k], $aux_start, $len;
		    $new_seqs[$k] = $str;
		}
		# Removes gap-only columns
		my $len = length($new_seqs[0]);
                my $size = scalar(@new_seqs);
                for(my $k = 0; $k < $len; ++$k){
                    my @col = @{get_column($k, \@new_seqs)};
                    my $gap = 100*(gap_percentage(\@col));
                    if($gap == 100){
                        for(my $l = 0; $l < $size; ++$l){
                            my @aux = split("", $new_seqs[$l]);
                            splice @aux,$k, 1;
                            $new_seqs[$l] = join("", @aux);
                        }
                        --$k;
                        --$len;
                    }
                }
                if($len < $min_size_block){
		    print LOG "Selected region is shorter than the minimum block size after removing gap-only columns: discarded\n";
                    next;
                }
  		my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@names, \@new_seqs, undef);
                @names = @{$aux_name};
                @new_seqs = @{$aux_seqs};
		my $nseq = scalar(@names); 
        	my $min_nseq = 2;
        	if(defined $amount){
            	    $min_nseq = $amount;
        	}
        	if($nseq < $min_nseq){
            	    print LOG "Number of sequences lower than $min_nseq: discarded\n";
            	    next;
        	}

		if(defined $gap_perc){
		    my @aux_seqs = ();
		    my @aux_names = ();
                    $len = ($end-$start)+1;
                    for(my $k = 0; $k < scalar(@names); ++$k){
                        my $j;
                        my @aux = split("", $new_seqs[$k]);
                        my $count = 0; 
			my $count_x = 0;
                        for($j = 0; $j < $len; ++$j){
                            if($aux[$j] eq '-'){
                                ++$count;
                            }
			    elsif(lc($aux[$j]) eq 'x'){
                                ++$count_x;
                            }
                        }
                        my $perc = ($count/$len)*100;
			my $perc_x = ($count_x/$len)*100;
                        if($perc <= $gap_perc and $perc_x <= $gap_perc){
			    push @aux_seqs, $new_seqs[$k];
			    push @aux_names, $names[$k];
                        }
			elsif($perc > $gap_perc and $perc_x > $gap_perc){
                            my $auxp = sprintf("%.2f", $perc);
                            my $auxx = sprintf("%.2f", $perc_x);
                            print LOG "Sequence $names[$k] removed ($auxp\% of gaps in sequence and $auxx\% of X\'s in sequence)\n";
                            next;
                        }
                        elsif($perc > $gap_perc){
                            my $auxp = sprintf("%.2f", $perc);
                            print LOG "Sequence $names[$k] removed ($auxp\% of gaps in sequence)\n";
                            next;
                        }
                        elsif($perc_x > $gap_perc){
                            my $auxx = sprintf("%.2f", $perc_x);
                            print LOG "Sequence $names[$k] removed ($auxx\% of X\'s in sequence)\n";
                            next;
                        }
                    }
		    my $nseq = scalar(@aux_names); #`grep -c ">" $seq_without_red`;
            	    my $min_nseq = 2;
            	    if(defined $amount){
                	$min_nseq = $amount;
            	    }
            	    if($nseq < $min_nseq){
                	print LOG "Number of sequences lower than $min_nseq: discarded \n";
                	next;
            	    }
		    else{
			open(FILE, ">$file_new");
			for(my $k = 0; $k < scalar(@aux_names); ++$k){
                            print FILE ">$aux_names[$k]\n";
                            print FILE "$aux_seqs[$k]\n";
                    	}
		     	close(FILE);
	 	    }
                }
                else{
		    open(FILE, ">$file_new");
                    for(my $k = 0; $k < scalar(@names); ++$k){
                        print FILE ">$names[$k]\n";
                        print FILE "$new_seqs[$k]\n";
                    }
		    close(FILE);
                }
            }
            else{
                print LOG "Selected region is shorter than the minimum block size: discarded\n";
            }
	    print LOG "\n";
	}	
    }
    print LOG "\n";
}
##########################################################################################
#                               Generates the final files                                #
##########################################################################################

sub generateFinalFiles{
    my $aux_blocks = shift;
    my $aux_values = shift;
    my $aux_scores = shift;
    my $dir = shift;
    my $category = shift;

    my @blocks = @{$aux_blocks};
    my @values = @{$aux_values};
    my @scores = @{$aux_scores};
    my @valid_blocks = ();
    for(my $i = 0; $i < scalar(@values); ++$i){
    	$valid_blocks[$i] = 0;
    }
    
    my $file_blocks;
    my $file;
    if(defined $category){
    	$file_blocks = $dir."/"."selected_blocks_".$category.".txt";
	$file = $dir."/scores_".$category.".csv";
    }
    else{
	$file = $dir."/scores.csv";
	$file_blocks = $dir."/"."selected_blocks.txt";
    }
    open(FILE, ">$file_blocks");
    print LOG "\nFinal valid blocks:\n";

    for(my $i = 0; $i < scalar(@blocks); ++$i){
    	my @aux = split("-", $blocks[$i]);
    	my $start = $aux[0];
    	my $end = $aux[1];
    	print LOG "$start-$end\n";
    	print FILE "$start-$end\n";
    	for(my $j = ($aux[0]-1); $j <= ($aux[1]-1); ++$j){
            $valid_blocks[$j] = $scores[$j];
    	}
    }
    close(FILE);
    
    printFileOnlyOption($file, \@scores, \@valid_blocks);

}

##########################################################################################
#                               Generates the final files                                #
##########################################################################################

sub generateFinalFilesMix{
    my $aux_blocks = shift;
    my $aux_blocks_mi = shift;
    my $aux_blocks_sh = shift;
    my $aux_values = shift;
    my $aux_values_mi = shift;
    my $aux_values_sh = shift;
    my $aux_scores = shift;
    my $dir = shift;
    my $category = shift;

    my @blocks = @{$aux_blocks};
    my @blocks_mi = @{$aux_blocks_mi};
    my @blocks_sh = @{$aux_blocks_sh};
    my @values = @{$aux_values};
    my @values_mi = @{$aux_values_mi};
    my @values_sh = @{$aux_values_sh};
    my @scores = @{$aux_scores};
    my @valid_blocks = ();
    for(my $i = 0; $i < scalar(@values); ++$i){
        $valid_blocks[$i] = 0;
    }

    my $file_blocks;
    my $file;
    if(defined $category){
	$file_blocks = $dir."/"."selected_blocks_".$category.".txt";
	$file = $dir."/scores_".$category.".csv";
    }
    else{
	$file_blocks = $dir."/"."selected_blocks.txt";
	$file = $dir."/scores.csv";
    }
    open(FILE, ">$file_blocks");
    print LOG "\nSelected blocks file: $file_blocks\n";
    print LOG "\nFinal valid blocks:\n";

    for(my $i = 0; $i < scalar(@blocks); ++$i){
        my @aux = split("-", $blocks[$i]);
        my $start = $aux[0];
        my $end = $aux[1];
        print LOG "$start-$end\n";
        print FILE "$start-$end\n";
        for(my $j = ($aux[0]-1); $j <= ($aux[1]-1); ++$j){
	    if($scores[$j] == -1000){
            	$valid_blocks[$j] = 0;
	    }
	    else{
		$valid_blocks[$j] = $scores[$j];
	    }
        }
    }
    close(FILE);

    my @mi_v = ();
    my @sh_v = ();
    for(my $i = 0; $i < scalar(@values); ++$i){
        $mi_v[$i] = 0;
        $sh_v[$i] = 0;
    }
    my @mi = @{setsIntersection(\@blocks_mi, \@blocks)};
    for(my $i = 0; $i < scalar(@mi); ++$i){
        my @aux = split("-", $mi[$i]);
        for(my $j = $aux[0]; $j <= $aux[1]; ++$j){
	    if($values_mi[$j] == -1000){
                $mi_v[$j] = 0;
            }
	    else{
            	$mi_v[$j] = $values_mi[$j];
	    }
        }
    }

    my @sh = @{setsIntersection(\@blocks_sh, \@blocks)};
    for(my $i = 0; $i < scalar(@sh); ++$i){
        my @aux = split("-", $sh[$i]);
        for(my $j = $aux[0]; $j <= $aux[1]; ++$j){
	    if($values_sh[$j] == -1000){
                $sh_v[$j] = 0;
            }
            else{
            	$sh_v[$j] = $values_sh[$j];
	    }
        }
    } 
    printFileMOption($file, \@scores, \@mi_v, \@sh_v, \@valid_blocks);
}

##########################################################################################
#                 Selects the best sub-region from an alignment region                   #
##########################################################################################

sub selectBestRegion{
    my $start = shift;
    my $end = shift;
    my $aux_values = shift;
    my $count = shift;
    my $aux_gaps = shift;
    my @value = @{$aux_values};    
    my @gaps_only = @{$aux_gaps};
    my @blocks = ();
    my $size = $end - $start + 1;
    for(my $i = 0; $i < $size; ++$i){
	my $s = $i;
      	my $e = ($i + $max_block_size+$count)-1;
      	if($e >= $size){
            last;
        }
      	my $sum = 0;
      	my $aux_start = $start + $i;
      	my $aux_end = $start + $e;
	my $win = 0;
	if($gaps_only[$aux_start] == 1){
	    next;
	}
      	for(my $j = $aux_start; $j < $aux_end; ++$j){	    
	    if($win >= $max_block_size){
		last;
	    }
	    if($aux_end >= $end){
		last;
	    }
	    if($gaps_only[$j] == 1){
		++$aux_end;
	    }	    
      	    elsif($value[$j] > 0){
            	$sum += $value[$j];
		++$win;
            }
      	}
	my $string = "$sum-$s-$aux_end";
        my $as = $aux_start + 1;
        my $ae = $aux_end + 1;
        print LOG "Region: $as - $ae: score sum=$sum\n";
        push @blocks, $string;
    }
    my @aux2 = split("-", $blocks[0]);
    my $sum_max = $aux2[0];
    my $start_max = $aux2[1];
    my $end_max = $aux2[2];
    my $t = scalar(@blocks);
    if(scalar(@blocks) > 0){
    	for(my $i = 1; $i < scalar(@blocks); ++$i){
    	    @aux2 = split("-", $blocks[$i]);
            if($aux2[0] > $sum_max){
            	$sum_max = $aux2[0];
            	$start_max = $aux2[1];
            	$end_max = $aux2[2];
            }
    	}
        my $as = $start_max + $start+ 1;
        my $ae = $end_max + 1; #$start + 1;
        print LOG "Top score region: $as - $ae: score sum=$sum_max\n";
        return ($as, $ae);
    }
    else{
	return(undef,undef);
    }
}

##########################################################################################
#                        Selects sequences from a alignment                              #
##########################################################################################

sub selectSequences{
    my $aux_aln = shift;
    my $start = shift;
    my $end = shift;
    my $a = shift;
    my $e = shift;
    my @new_seqs = ();
    my @alns = @{$aux_aln};
    for(my $k = 0; $k < scalar(@alns); ++$k){
    	my $start_a = ($a - 1) - $start;
        my $end_a = $end - ($e - 1);	
        my $len = ($e - 1) - ($a - 1) + 1;
	if($start_a < 0){
	    $start_a = 0;
	}
        my $str = substr $alns[$k], $start_a, $len;
        $new_seqs[$k] = $str;
    }
    my $len = length($new_seqs[0]);
    my $size = scalar(@new_seqs);
    my @del = ();
    for(my $k = 0; $k < $len; ++$k){
     	my @col = @{get_column($k, \@new_seqs)};
        my $gap = 100*(gap_percentage(\@col));
        if($gap == 100){
	    push @del, $k;
	}
    }
    for(my $k = 0; $k < $size; ++$k){
	my @aux = split("", $new_seqs[$k]);
	my @aux2 = ();
   	for(my $d = 0; $d < scalar(@aux); ++$d){
	    if(contains($d, \@del) == 0){
		push @aux2, $aux[$d];
	    }
	}
	$new_seqs[$k] = join("", @aux2);
    } 

    return \@new_seqs; 
}

##########################################################################################
#                         Verifies if a vector contains a value                          #
##########################################################################################

sub contains {
    my $num = shift;
    my $aux = shift;
    my @vector = @{$aux};
    for(my $i = 0; $i < scalar(@vector); ++$i){
	if($vector[$i] == $num){
	    return 1;	
	}
    }
    return 0;
}

##########################################################################################
#                                Calculates the gap percentage                           #
##########################################################################################

sub calculateGapPercentage{
    my $aux_names = shift;
    my $aux_alns = shift;
    my $start = shift;
    my $end = shift;
    my $len = ($end-$start)+1;
    my @aux_n = ();
    my @aux_s = ();
    my @nms = @{$aux_names};
    my @alns = @{$aux_alns};
    for(my $k = 0; $k < scalar(@nms); ++$k){
        my $j;
        my @aux = split("", $alns[$k]);
        my $count = 0;
        my $count_x = 0;
        for($j = 0; $j < $len; ++$j){
            if($aux[$j] eq '-'){
            	++$count;
            }
            elsif(lc($aux[$j]) eq 'x'){
                ++$count_x
            }
        }
	my $perc;
        my $perc_x;
	if($len  == 0){
	    $perc = 0;
	    $perc_x = 0;
  	}
	else{
	    $perc = ($count/$len)*100;
	    $perc_x = ($count_x/$len)*100;
	}
        if($perc <= $gap_perc and $perc_x <= $gap_perc){
            push @aux_n, $nms[$k];
            push @aux_s, $alns[$k];
        }
        elsif($perc > $gap_perc and $perc_x > $gap_perc){
             my $auxp = sprintf("%.2f", $perc);
             my $auxx = sprintf("%.2f", $perc_x);
             print LOG "Sequence $nms[$k] removed ($auxp\% of gaps in sequence and $auxx\% of X\'s in sequence)\n";
             next;
        }
        elsif($perc > $gap_perc){
	     my $auxp = sprintf("%.2f", $perc);
             print LOG "Sequence $nms[$k] removed ($auxp\% of gaps in sequence)\n";
             next;
        }
        elsif($perc_x > $gap_perc){
             my $auxx = sprintf("%.2f", $perc_x);
             print LOG "Sequence $nms[$k] removed ($auxx\% of X\'s in sequence)\n";
             next;
        }
    }
    return \@aux_n, \@aux_s;
}

##########################################################################################
#                Searchs for more sub-regions in a selected alignment block              #
##########################################################################################

sub searchForMoreRegions{
    my $original_start = shift;
    my $original_end = shift;
    my $block_start = shift;
    my $block_end = shift;
    my $aux_values = shift;
    my $aux_gaps = shift;
    my @values = @{$aux_values};
    my @gaps_only = @{$aux_gaps};
    my $str = searchRecursively ($original_start, $original_end, $block_start, $block_end, \@values, \@gaps_only);
    my @aux = split("_", $str);
    @aux = @{sortCoordinates(\@aux)};
    my $index;
    for(my $i = 0; $i < scalar(@aux); ++$i){
	my @aux_cood = split("-", $aux[$i]);
	if($aux_cood[0] == $block_start and $aux_cood[1] == $block_end){
	    $index = $i;
	    last;
	}
    }
    splice @aux, $index, 1;
    return \@aux;
}

##########################################################################################
#                              Searches recursively for regions                          #
##########################################################################################

sub searchRecursively{
    my $original_start = shift;
    my $original_end = shift;
    my $block_start = shift;
    my $block_end = shift;
    my $aux_values = shift;
    my $aux_gaps = shift;
    my $aux_coord = $original_start+1;
    my @values = @{$aux_values};
    my @gaps_only = @{$aux_gaps};
    my $dist_left = $block_start - $original_start;
    my $dist_right = $original_end - $block_end;
    if($dist_left < $min_size_block and $dist_right < $min_size_block){ 
	return "$block_start-$block_end";
    }
    elsif($dist_left >= $max_block_size and $dist_right >= $max_block_size){	
  	my ($start_left, $end_left) = selectBestRegion($original_start, $block_start, \@values, 0, \@gaps_only);
	my ($start_right, $end_right) = selectBestRegion($block_end, $original_end, \@values, 0, \@gaps_only);
	if((defined $start_left) and (defined $end_left) and (defined $start_right) and (defined $end_right)){
            return searchRecursively($original_start, $block_start, $start_left, $end_left, \@values, \@gaps_only)."_"."$block_start-$block_end"."_".searchRecursively($block_end, $original_end, $start_right, $end_right, \@values, \@gaps_only);
	}
	elsif((defined $start_left) and (defined $end_left)){
	    return searchRecursively($original_start, $block_start, $start_left, $end_left, \@values, \@gaps_only)."_"."$block_start-$block_end";
	}
	elsif((defined $start_right) and (defined $end_right)){
	    return "$block_start-$block_end"."_".searchRecursively($block_end, $original_end, $start_right, $end_right, \@values, \@gaps_only);
	}
	else{
	    return "$block_start-$block_end";
	}
    }
    elsif(($dist_left < $max_block_size and $dist_left >= $min_size_block) and ($dist_right < $max_block_size and $dist_right >= $min_size_block)){
	return "$block_start-$block_end"."_"."$original_start-$block_start"."_"."$block_end-$original_end";
    }
    elsif($dist_left < $min_size_block and $dist_right > $max_block_size){
	my ($start_right, $end_right) = selectBestRegion($block_end, $original_end, \@values, 0, \@gaps_only);
	if((defined $start_right) and (defined $end_right)){
	    return "$block_start-$block_end"."_".searchRecursively($block_end, $original_end, $start_right, $end_right, \@values, \@gaps_only);
	}
	else{
	    return "$block_start-$block_end";
	}
    }
    elsif($dist_left < $min_size_block and ($dist_right <= $max_block_size and $dist_right >= $min_size_block)){
        return "$block_start-$block_end"."_"."$block_end-$original_end";
    }    
    elsif($dist_right <  $min_size_block and $dist_left > $max_block_size){
        my ($start_left, $end_left) = selectBestRegion($original_start, $block_start, \@values, 0, \@gaps_only);
	if((defined $start_left) and (defined $end_left)){
            return "$block_start-$block_end"."_".searchRecursively($original_start, $block_start, $start_left, $end_left, \@values, \@gaps_only);
	}
	else{
	    return "$block_start-$block_end";
	}
    }
    elsif($dist_right < $min_size_block and ($dist_left <= $max_block_size and $dist_left >= $min_size_block)){
        return "$block_start-$block_end"."_"."$original_start-$block_start";
    }
    elsif(($dist_left <= $max_block_size and $dist_left >= $min_size_block) and ($dist_right > $max_block_size)){
        my ($start_right, $end_right) = selectBestRegion($block_end, $original_end, \@values, 0, \@gaps_only);	
	if((defined $start_right) and (defined $end_right)){
            return "$block_start-$block_end"."_"."$original_start-$block_start"."_".searchRecursively($block_end, $original_end, $start_right, $end_right, \@values, \@gaps_only);
	}
	else{
	    return "$block_start-$block_end";
	}
    }
    elsif(($dist_right <= $max_block_size and $dist_right >= $min_size_block) and ($dist_left > $max_block_size)){
	my ($start_left, $end_left) = selectBestRegion($original_start, $block_start, \@values, 0, \@gaps_only);
	if((defined $start_left) and (defined $end_left)){
            return "$block_start-$block_end"."_".searchRecursively($original_start, $block_start, $start_left, $end_left, \@values, \@gaps_only)."_"."$block_end-$original_end";
	}
	else{
	    return "$block_end-$original_end";
	}
    }
}

##########################################################################################
#                                  Sort block coordinates                               #
##########################################################################################

sub sortCoordinates{
    my $aux_vector = shift;
    my @vector = @{$aux_vector};
    my $len = scalar(@vector);
    for(my $i = 0; $i < $len - 1; ++$i){
        my @aux = split("-", $vector[$i]);
        my $start = $aux[0];
        my $index = $i;
        my $lower = $start;
        for(my $j = $i+1; $j < $len; ++$j){
            my @aux2 = split("-", $vector[$j]);
            my $start2 = $aux2[0];
            if($lower > $start2){
                $lower = $start2;
                $index = $j;
            }
        }
        if($index != $i){
            my $aux3 = $vector[$i];
            $vector[$i] = $vector[$index];
            $vector[$index] = $aux3;
        }
    }
    return \@vector;
}

##########################################################################################
#                         Generates a FASTA file from a set of sequences                 #
##########################################################################################

sub generateRegion{
    my $original_start = shift;
    my $original_end = shift;
    my $block_start = shift;
    my $block_end = shift;
    my $aux_value = shift;
    my $aux_nms = shift;
    my $aux_alns = shift;
    my $dir = shift;
    my $num_seq = shift;
    my $category = shift;
    my @values = @{$aux_value};
    my @names = @{$aux_nms};
    my @alns = @{$aux_alns};
    my $file_new;
    my $a = $block_start;
    my $e = $block_end;
    if($measure eq "c"){
    	if($input_file_prefix =~ /_/){
            my @aux = split("_", $input_file_prefix);
            $file_new =  $dir."/".$aux[0]."_".$num_seq."_".$a."-".$e.".fasta";
        }
        else{
            $file_new =  $dir."/".$category."_".$num_seq."_".$a."-".$e.".fasta";
        }
    }
    else{
        $file_new =  $dir."/".$category."_".$num_seq."_".$a."-".$e.".fasta";
    }
    my @new_seqs = @{selectSequences(\@alns, $original_start, $original_end, $block_start, $block_end)};    
    my $len = length($new_seqs[0]);
    if($len < $min_size_block){
    	print LOG "Selected region is shorter than the minimum block size after removing gap-only columns: discarded\n\n";
        return -1;
    }
    my ($count, $aux_name, $aux_seqs) = remove_redundancy(\@names, \@new_seqs, undef);
    @names = @{$aux_name};
    @new_seqs = @{$aux_seqs};
    my $nseq = scalar(@names);
    my $min_nseq = 2;
    if(defined $amount){
    	$min_nseq = $amount;
    }
    if($nseq < $min_nseq){
    	print LOG "Number of sequences lower than $min_nseq: discarded\n\n";
        return -1;
     }
     if(defined $gap_perc){
	my $start_a = $block_start - $original_start;
        my $end_a = $original_end - $block_end;
     	($aux_name, $aux_seqs) = calculateGapPercentage(\@names, \@new_seqs, $start_a, $end_a);
        my @aux_n = @{$aux_name};
        my @aux_s = @{$aux_seqs};
        my $min_nseq = 2;
        if(defined $amount){
            $min_nseq = $amount;
        }
        if(scalar(@aux_n) >= $min_nseq){
            open(FILE, ">$file_new");
            for(my $k = 0; $k < scalar(@aux_n); ++$k){
            	print FILE ">$aux_n[$k]\n";
                print FILE "$aux_s[$k]\n";
            }
            close(FILE);
         }
    }
    else{
        open(FILE, ">$file_new");
	for(my $k = 0; $k < scalar(@names); ++$k){
            print FILE ">$names[$k]\n";
            print FILE "$new_seqs[$k]\n";
        }
        close(FILE);
    }
    return 0; 	
}

##########################################################################################
#                               Selects gap-only columns                                 #
##########################################################################################

sub get_gapOnly{
   my $aux = shift;
   my @seqs = @{$aux};
   my $tam = scalar(@seqs);
   my @matrix = ();
   my $i, $j;
   my $num_col = length($seqs[0]);
   my @gaps_only = ();
   for($i = 0; $i < $tam; ++$i){
        my @aux = split('', $seqs[$i]);
        for($j = 0; $j < $num_col; ++$j){
            $matrix[$i][$j] = $aux[$j];
        }
    }
    for($i = 0; $i < $num_col; ++$i){
        my $count = 1;
        for($j = 0; $j < $tam; ++$j){
            if($matrix[$j][$i] ne "-"){
                $count = 0;
                last;
            }
        }
        $gaps_only[$i] = $count;
    }
    return \@gaps_only;
}
