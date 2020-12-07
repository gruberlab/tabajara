# TABAJARA - Tool for Alignment Block Analysis Joining Appropriate Rational Approaches

TABAJARA is a computational tool for rational design of profile HMMs. Starting from a multiple sequence alignment (MSA), TABAJARA is able to find blocks that are either conserved across all sequences or discriminative for two specific groups of sequences. 

#   Instalation

Tabajara does not need to be installed. The user should only download the tabajara.pl file.

# Requirements

TABAJARA requires the program hmmbuild (HMMER3 package - (http://hmmer.org/) to build profile HMMs. The program must be located in a directory listed in the PATH of the operating system.

# Usage
```
perl tabajara.pl -i <input_file> -o <output_directory>
```
### Mandatory parameters:
```
-b|minimum_block_size <integer>         Minimum block size (must be ≥ w).
  
-i|input_file <file name>               Input file (multiple alignment in FASTA or CLustal formats).
  
-fl|full_length <yes|no>                Use full-length sequence for model construction (default = no).

##### IMPORTANT: When using -fl no (default), the following parameters are also mandatory:

-m|method <c|d>                         Method to select blocks. Options:
                                         c - Conservation
                                         d - Discrimination

-t|threshold <decimal>                  Threshold of the alignment position for block extraction  (valid values: 0 to 1).

-p|percentage <integer>                 Percentage of positions in sliding window with score ≥ t.

-w|window_size <integer>                Window size for block extraction.
```

### Optional parameters:
```
-c|category <string>                    Name of major category to be analysed.

-conf                                   Configuration file

-clean <yes|no>                         Remove all intermediate files during the execution of Tabajara (default = yes).

-cs|cutoff_score <yes|no>               Insert cutoff scores in the profile HMMs (default = no).

-di|discard_sequences <yes|no>          Discard identical sequences (default = no).

-gc|gap_cutoff <integer>                Gap cutoff. Do not score columns that contain more than gap cutoff fraction gaps (default = 30).

-gs|gap_sequence <decimal>              Percentage of allowed gaps in each sequence.

-hb|hmmbuild_parameters <string>        Any set of hmmbuild's valid parameters can be entered under double quotes.
                                            Example: -hb \"--wblosum --wid 0.8\". Default: no parameters.
-h|help                                 Help.

-mb|maximum_block_size <integer>        Maximum block size.

-md|maximum_distance <integer>          Maximum accepted distance for two alignment blocks to be joined as a single block.

-o|output                               Output directory.

-pc|minimum_category_rate               Minimum percentage of categories in the training set that must meet the criteria defined by the parameter pd for profile HMMs built in conservation mode (default = N).

-pd|minimum_detection_rate <decimal>    Minimum detection rate (in percentage) of training set (MSA) sequences by the constructed HMM.
                                        This is equivalent to the minimum accepted sensitivity of the model (default = 80).

-pt|maximum_percentage_ratio <decimal>  This parameter applies for the validation of HMMs designed for group discrimination. Given two groups,
                                        A (selected category) and B (remaining sequences), an hmmsearch is performed with the constructed model
                                        against the training set (MSA) sequences. The highest score obtained by group B sequences must be lower
                                        than n % (defined by parameter -pt) of the lowest score obtained by group A sequences, otherwise the model
                                        is not accepted (default = 80).

-sa|sequence_amount <integer>           Minimum sequence amount of MSA to construct hmm (default = 2).

-sr|search_more_regions <yes|no>        This parameter is used to search for more subregions in a previously selected region (default = no).

-sv|score_value <decimal>               This parameter is used to calculate the minimum score for HMMs designed by Tabajara. The minimum score will be used as                                               suggested score if calculated score is lower than this score. This minimum score is calculated by
                                        mutiplying the value of -sv by HMM length (default = 1).

-v|version                              Version.

-wg|discard_windows <yes|no>            Discard sliding windows presenting gaps (default = no).
``` 
# HMM-prospector: A tool for surveying genomic and metagenomic data with profile HMMs

HMM-Prospector is a Perl script that uses a single or multiple profile HMM file as a query in similarity searches against a FASTQ/FASTA dataset using hmmsearch program. HMM-Prospector processes the results and generates tabular files with qualitative and quantitative results. Previously run hmmsearch result files (short tabular format) can also be used as datasets.

#   Instalation

HMM-Prospector also does not need to be installed. The user should only download the hmm-prospector.pl file.

# Requirements

- hmmsearch: HMMER3package-http://hmmer.org/. 

- transeq: EMBOSS package - http://emboss.sourceforge.net/ .

- fastq_to_fasta: FASTX-Toolkit–http://hannonlab.cshl.edu/fastx_toolkit/.

P.S: All programs must be located in a directory listed in the PATH of the operating system.

# Usage

```
perl hmm-prospector.pl -d <file> -s|-e <decimal>  <optional parameters>
```  

### Mandatory parameters:

```
-d  <file>        : Dataset (FASTQ, FASTA or hmmsearch's tabular output file)
```

### OPTIONAL PARAMETERS:

```
-a                : Directory containing profile HMM annotations (valid only when using vFam models as input).

-cpu              : Number of threads to be used by hmmsearch.

-e | -s <decimal> : E-value (-e) or score (-s) threshold value. Report hmmsearch hits that present values equal to
		    or lower than E-value or equal to or larger than score. One parameter and the respective value
		    must be provided. If an hmmsearch tabular result file is used as input, then parameters -e or -s
		    become mandatory.

-h|help           : Show this help message.

-i                : Input file (single or multiple profile HMMs) - mandatory when using a FASTA or FASTQ dataset.

-o             	  : Output directory (default = output_dir).

-r 		  : Ignore cutoff scores in the profile HMMs and use a custom value defined by parameters -e or -s
		    for all input models (default = yes). If -r no is used, HMM-Prospector will use the cutoff scores
		    specified in the respective CUTOFF SCORE tag of each profile HMM. For models not containing cutoff
		    values, HMM-Prospector will use the cutoff value specified by the parameter -e ou -s. If none of these
		    parameters has been specified, the program will then use hmmsearch's default cutoff value (-E 10).

-v|version        : Version.
```

# Contact

To report bugs, to ask for help and to give any feedback, please contact Arthur Gruber (argruber@usp.br) or Liliane S. Oliveira (liliane.sntn@gmail.com).
