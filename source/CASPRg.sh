#!/bin/bash

###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: generate_counts.sh [-h or --help]
                          [-f or --fastq-forward]
                          [-r or --fastq-reverse]
                          [-l or --library]
                          [-e or --exper-design]
                          [-c or --controls]
                          [-o or --orientation]
                          [-a or --adapter-f]
                          [-A or --adapter-r]
                          [-m or --mismatches]
                          [-b or --bases-aligned]
                          [-y or --fdr-threshold]
                          [-t or --threads]
                          [-k or --keep-tmp]
                          [-q or --output-folder]
"""
}
function print_error {
  printf "Error with the parameters"
  print_usage
  echo "Try -h or --help for more information"
}
function print_help {
  print_usage
  printf """
Optional arguments:
    -h, --help:
                Show this help message and exit.
    -r, --fastq-reverse:
                fastq files with reverse reads (only for pgRNAs).
                Accepted formats: fastq or fastq.gz.
                Required name: SampleName_reverse.fastq(.gz).
                Use commas to separate files; write them in quotes.
                Example: \"ctrl1_reverse.fastq, d1_reverse.fastq\".
                If not provided, only sgRNAs will be considered.
    -c, --controls:
                txt file with 2 columns.
                1st column: name of the control genes.
                2nd column: type of control.
                Use: Positive, Negative or Neutral.
    -o, --orientation:
                Orientation of the paired guides in plasmid vector.
                Use it only for pgRNAs.
                If 5'-gRNAs-3': write 53. If 3'-gRNAs-5': write 35.
                Default: 53.
    -a, --adapter-f:
                Sequence of the adapter in the forward read.
                It must contain a sequence of nucleotides.
                Default: ACCG.
    -A, --adapter-r:
                Sequence of the adapter in the reverse read.
                Use it only for pgRNAs.
                It must contain a sequence of nucleotides.
                Default: AAAC.
    -m, --mismatches:
                Maximum number of mismatches tolerated for counts.
                It must be an integer.
                Default: 0.
    -b, --bases-aligned:
                Minimum number of bases that must be aligned for counts.
                It must be an integer.
                Default: 20 for sgRNAs and 35 for pgRNAs.
    -y, --fdr-threshold:
                FDR threshold for the identification of hits.
                Must be a number between 0 and 1 (limits not included).
                Default: 0.1.
    -t, --threads:
                Number of threads that will be used.
                It must be an integer.
                Set it to the number of available cores.
                Default: 8.
    -k, --keep-tmp:
                If used, keep the intermediate files.
                Default: remove them.
    -q, --output-dir:
                Path and name of the directory to save outputs.
                The directory must exist. It will not be created.
                Default: current directory.

Required arguments:
    -f, --fastq-forward:
                fastq files with forward reads.
                Accepted formats: fastq or fastq.gz.
                Required name: SampleName_forward.fastq(.gz).
                Use commas to separate files; write them in quotes.
                Example: \"ctrl1_forward.fastq, d1_forward.fastq\".
    -l, --library:
                txt file with 2 or 3 columns, without header.
                1st column: IDs of the guide RNAs. Required format:
                NameOfGene_gRNAtag. Example: BRCA1_guide1.
                2nd column: Sequence of the 1st guide RNA.
                3nd column: Sequence of the 2nd guide RNA, if pgRNAs.
    -e, --exper-design:
                txt file with 3 columns, without header.
                1st column: Eeach row contains a SampleName.
                SampleName's must correspond to the fastq files names.
                2nd column: Index (integer) showing which samples
                (replicates) are the same in different periods of time.
                3rd column: Condition of each sample. Required:
                \"control\" or \"treated\" followed by a number.
                Example:    ctr1  1  control1
                            ctr2  2  control1
                            d1    1  treated1
                            d2    2  treated1
                            ctr1  1  control2
                            ctr3  3  control2
                            d1    1  treated2
                            d3    3  treated2
                --> Analysis 1: ctr1 & ctr2 vs d1 & d2
                --> Analysis 2: ctr1 & ctr3 vs d1 $ d3

Important:
    The number of input files in -f and -r must be the same (for pgRNAs).
    The SampleName's of -f and -r files must be identical.
    Forward and reverse fastq files must be provided in the same order.
"""
}

# Define inputs
for ARGS in "$@"; do
  shift
        case "$ARGS" in
                "--fastq-reverse") set -- "$@" "-r" ;;
                "--fastq-forward") set -- "$@" "-f" ;;
                "--library") set -- "$@" "-l" ;;
                "--exper-design") set -- "$@" "-e" ;;
                "--controls") set -- "$@" "-c" ;;
                "--orientation") set -- "$@" "-o" ;;
                "--adapter-f") set -- "$@" "-a" ;;
                "--adapter-r") set -- "$@" "-A" ;;
                "--mismatches") set -- "$@" "-m" ;;
                "--bases-aligned") set -- "$@" "-b" ;;
                "--fdr-threshold") set -- "$@" "-y" ;;
                "--threads") set -- "$@" "-t" ;;
                "--keep-tmp") set -- "$@" "-k" ;;
                "--output-dir") set -- "$@" "-q" ;;
                "--help") set -- "$@" "-h" ;;
                *) set - "$@" "$ARGS"
        esac
done

# Define defaults
o=53; a="ACCG"; A="AAAC"; m=0; b=0; t=1; k=0; y=0.1; q="./"

# Define all parameters
while getopts 'r:f:l:e:c:o::a::A::m::b::y::t::q::kh' flag; do
        case "${flag}" in
                r) r=${OPTARG} ;;
                f) f=${OPTARG} ;;
                l) l=${OPTARG} ;;
                e) e=${OPTARG} ;;
                c) c=${OPTARG} ;;
                o) o=${OPTARG} ;;
                a) a=${OPTARG} ;;
                A) A=${OPTARG} ;;
                m) m=${OPTARG} ;;
                b) b=${OPTARG} ;;
                y) y=${OPTARG} ;;
                t) t=${OPTARG} ;;
                q) q=${OPTARG} ;;
                k) k=1 ;;
                h) print_help
                   exit 1;;
                *) print_error
                    exit 1;;
        esac
done

# Modify inputs to work with them in an easier way
r=$(echo $r | tr "," " ")
f=$(echo $f | tr "," " ")

# If single guides are used, set --orientation always to 53
if [[ $r == "" ]]; then
  o=53
  if [[ $b == 0 ]]; then b=20; fi
else
  if [[ $b == 0 ]]; then b=35; fi
fi

########################################
## Create folders to save final files ##
########################################

# Check that output directory exists
ls $q >/dev/null 2>&1
if [[ $(echo $?) != 0 ]]; then # Exit if there has been an error
  echo "Error: output directory doesn't exist"
  exit 1
fi
# Make sure that last character of output directory is "/"
if [[ $(echo "${q: -1}") != "/" ]]; then
  q="$q/"
fi

rm -rf "${q}intermediate" && \
mkdir "${q}intermediate" && chmod +xwr "${q}intermediate"
rm -rf "${q}outputs" && mkdir "${q}outputs" && chmod +xwr "${q}outputs"

########################################
# Execute the needed scripts in order ##
########################################

# Get the directory of this file to be able to call the others
currentdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Display any error that is initially detected
bash $currentdir/basic_errors.sh \
"$f" "$r" $l $e $o $a $A $m $b $y $t $q "$c" # Inputs
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.

# Create a fake genome with the guide RNAs
bash $currentdir/fake_genome.sh $l $t $q 2> "${q}intermediate/err.txt"
if [[ $(cat "${q}intermediate/err.txt") != "" ]]; then # Exit if error
  echo "Error: problems with the library file. Check it again."
  echo "See ${q}intermediate/err.txt for more information"
  exit 1
fi

# Perform QC of the fastq files. No new errors are expected in this step
bash $currentdir/quality_control.sh \
      "$f" "$r" $t $q 2> "${q}intermediate/err.txt"

# Trimming and arrangement of the reads
bash $currentdir/trimming_reads.sh \
          "$f" "$r" $l $o $a $A $q $t 2> "${q}intermediate/err.txt"
if [[ $(cat "${q}intermediate/err.txt") != "" ]]; then # Exit if error
  echo "Error: problems with the input parameters. Check them again."
  echo "See ${q}intermediate/err.txt for more information"
  exit 1
fi

# Create plots with the trimming statistics
printf "\nGenerating plots to see trimming statistics\n"
# Collect needed information previously created
triminf=$(cat "${q}intermediate/useful_information.txt" | \
          awk 'NR>1' | sort -k1 | awk '{print $2}')
# Run R script
Rscript --vanilla $currentdir/trimming_statistics.R \
          "$triminf" "$q" ${q}intermediate/trim_stat* \
|| (echo "Problem with R. Check if the R version is correct." && exit 1)
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.
echo "Plots were created successfully"

# Align the reads to the fake genome and count the number of reads per guide
bash $currentdir/alignment_counts.sh \
      $m $b $t $q "$r" $l 2> "${q}intermediate/err.txt"
if [[ $(cat "${q}intermediate/err.txt") != "" ]]; then # Exit if error
  echo "Error: problems with the input parameters. Check them again."
  echo "See ${q}intermediate/err.txt for more information."
  exit 1
fi

# Create plots with the alignment statistics
printf "\nGenerating plots to see alignment statistics\n"
Rscript --vanilla $currentdir/alignment_statistics.R \
        $q ${q}intermediate/Statistics_alignment* \
|| (echo "Problem with R. Check if the R version is correct." && exit 1)
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.

# Create plots to see what happens with unmapped reads
Rscript --vanilla $currentdir/alignment_unmapped.R \
        $q ${q}intermediate/Statistics_unmapped* \
|| (echo "Problem with R. Check if the R version is correct." && exit 1)
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.
echo "Plots were created successfully"

# Compute the test to see the significant genes.
bash \
  $currentdir/test.sh $e $y $q $currentdir "$c" 2> "${q}intermediate/err.txt"
if [[ $(echo $?) == 2 ]]; then exit 1; fi # Exit if there has been an error.
if [[ $(cat ${q}intermediate/err.txt | grep -oEi "^ERROR") != "" ]]; then
  echo "Error: problems with MAGeCK."
  echo "See ${q}intermediate/err.txt for more information"
  exit 1
fi

# Visualize some of the results at the guides level
printf "\nGenerating plots to visualize general results\n"
Rscript --vanilla $currentdir/create_graphs.R \
        $q "$c" ${q}intermediate/*sgrna_summary.txt \
|| (echo "Problem with R. Check if the R version is correct." && exit 1)
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.
echo "Plots were created successfully"

# Create files for the visualization with VISPR
bash $currentdir/visualization.sh ${q} $currentdir

###############################
## Remove intermediate files ##
###############################

if [[ $k == 0 ]]; then rm -r ${q}intermediate; fi # Remove intermediate files

##########
## DONE ##
##########
