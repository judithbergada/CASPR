#!/bin/bash

#########################
## Define basic errors ##
#########################

# Relationship between input parameters and the ones used here
f=$1; r=$2; l=$3; e=$4; o=$5; a=$6; A=$7
m=$8; b=$9; y=${10}; t=${11}; q=${12}; c=${13}
start=${14}; pause=${15}

printf "\nChecking that the format of inputs is correct\n"

# Check if all of the required modules are installed and
# show information on how to install the missing modules.
if ! [ -x "$(command -v fastqc)" ]; then
  printf "Missing: "
  echo "fastqc not found"
  echo "Information on the installation:"
  echo "http://www.bioinformatics.babraham.ac.uk/projects/download.html"
  exit 1
fi
if ! [ -x "$(command -v cutadapt)" ]; then
  printf "Missing: "
  echo "cutadapt not found"
  echo "Information on the installation:"
  echo "https://cutadapt.readthedocs.io/en/stable/installation.html"
  exit 1
fi
if ! [ -x "$(command -v reformat.sh)" ]; then
  printf "Missing: "
  echo "bbmap not found"
  echo "Information on the installation:"
  echo "https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/"
  exit 1
fi
if ! [ -x "$(command -v STAR)" ]; then
  printf "Missing: "
  echo "STAR not found"
  echo "Information on the installation:"
  echo "https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
  exit 1
fi
if ! [ -x "$(command -v samtools)" ]; then
  printf "Missing: "
  echo "samtools not found"
  echo "Information on the installation:"
  echo "http://www.htslib.org/download/"
  exit 1
fi
if ! [ -x "$(command -v gs)" ]; then
  printf "Missing: "
  echo "Ghostscript not found"
  echo "Information on the installation:"
  echo "https://www.ghostscript.com/doc/9.20/Make.htm"
  exit 1
fi
if ! [ -x "$(command -v mageck)" ]; then
  printf "Missing: "
  echo "MAGeCK not found"
  echo "Information on the installation:"
  echo "https://bitbucket.org/liulab/mageck-vispr"
  exit 1
fi
if ! [ -x "$(command -v vispr)" ]; then
  printf "Missing: "
  echo "VISPR not found"
  echo "Information on the installation:"
  echo "https://bitbucket.org/liulab/mageck-vispr"
  exit 1
fi
if ! [ -x "$(command -v R)" ]; then
  printf "Missing: "
  echo "R not found"
  echo "Information on the installation:"
  echo "https://www.r-project.org/"
  exit 1
fi
echo "All of the required programs are properly installed"

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

###########################################################
## Check if the inputs are given in the expected formats ##
###########################################################

# Check that the forward files are indeed provided
if [[ $f == "" ]]; then
  echo "Error: forward fastq files are required"
  exit 1
fi

# Check if the library file is given as input and if it exists
if [ -z $l ]; then
  echo "Error: library is missing. Check the format of inputs"
  exit 1
fi
if [ ! -f $l ]; then
  echo "Error: the library file doesn't exist"
  exit 1
fi
# Check if the library file has hidden characters at the end, and remove them
cat $l | tr -d "\r" > ${q}/libnotused.txt
mv ${q}/libnotused.txt $l
# Replace any white space between columns by a tab
cat $l | tr "[:blank:]" "\t" > ${q}/libnotused.txt
mv ${q}/libnotused.txt $l
# Count the number of columns of the library
ncolslib=$(head -n1 $l | awk '{print NF}')

# Check if the experimental design file is given as input and if it exists
if [ -z $e ]; then
  echo "Error: experimental design is missing. Check the format of inputs"
  exit 1
fi
if [ ! -f $e ]; then
  echo "Error: the experimental design file doesn't exist"
  exit 1
fi
# Check if the exper.design file has hidden characters; remove them
cat $e | tr -d "\r" > ${q}/expnotused.txt
mv ${q}/expnotused.txt $e
# Replace any white space between columns by a tab
cat $e | tr "[:blank:]" "\t" > ${q}/expnotused.txt
mv ${q}/expnotused.txt $e
# Count the number of columns of the experimental design
ncolsexp=$(head -n1 $e | awk '{print NF}')
if [[ $ncolsexp != 3 ]]; then
  echo "Error: the experimental design file doesn't have 3 columns"
  exit 1
fi

# Check if the controls file exists (only if it is given as input)
if [[ $c != "" ]]; then
  if [ ! -f $c ]; then
    echo "Error: the controls file doesn't exist"
    exit 1
  fi
  # Check if the controls file has hidden characters; remove them
  cat $c | tr -d "\r" > ${q}/ctrlsnotused.txt
  mv ${q}/ctrlsnotused.txt $c
  # Replace any white space between columns by a tab
  cat $c | tr "[:blank:]" "\t" > ${q}/ctrlsnotused.txt
  mv ${q}/ctrlsnotused.txt $c
  # Count the number of columns of the controls file
  ncolscontrol=$(head -n1 $c | awk '{print NF}')
  if [[ $ncolscontrol != 2 ]]; then
    echo "Error: the controls file doesn't have 2 columns"
    exit 1
  fi
fi

# Check that -r fastq files and the number of columns in library are correct
if [[ $r != "" ]]; then
  for fastqfile in $r; do
    # Check that the fastq files do exist
    if [ ! -f ${fastqfile} ]; then
      echo "Error: some of the reverse fastq files don't exist"
      exit 1
    fi
    # Check that the library has 4 columns
    if [[ $ncolslib != 4 ]]; then
      echo "Error: the library doesn't have 4 columns, as needed for pgRNAs"
      exit 1
    fi
  done
  echo "The reverse fastq files are correct"
else
  if [[ $ncolslib != 3 ]]; then
    echo "Error: the library doesn't have 3 columns, as needed for sgRNAs"
    exit 1
  fi
fi

# Check that -f fastq files are correct
for fastqfile in $f; do
  # Check that the fastq files do exist
  if [ ! -f ${fastqfile} ]; then
    echo "Error: some of the forward fastq files don't exist"
    exit 1
  fi
done
echo "The forward fastq files are correct"

# Check that the number of forward and reverse fastq files is the same
if [[ $r != "" ]]; then
  rev=$(echo $r | wc -w)
  fwd=$(echo $f | wc -w)
  if [[ ! $rev -eq $fwd ]]; then
    echo "Error: the number of forward and reverse fastq files is different"
    exit 1
  fi
  echo "The number of forward and reverse fastq files is the same"
fi

# Check that -a, -A are a sequence of nucleotides and -m, -t are integers.
ad="^[ATCGatcg]+$"
val="^[0-9]+$"
if [[ ! $a =~ $ad || ! $A =~ $ad || \
      ! $m =~ $val || ! $b =~ $val || ! $t =~ $val ]]; then
  echo "Error: problem with the format of the inputs -a, -A, -m, -b or -t"
  exit 1
fi

# Check that -y is a number between zero and 1
limits="^0\.[0-9]+$"
if [[ ! $y =~ $limits ]]; then
  echo "Error: problem with the fdr-threshold. Must be between 0 and 1"
  echo "Remember: 0 and 1 are not included."
  exit 1
fi

# Check that -o is either 35 or 53.
if [[ $o != 35 && $o != 53 ]]; then
  echo "Error: problem with the format of the orientation. Must be 35 or 53"
  exit 1
fi

# Check that --start is one of the possible options
if [[ $start != "" && $start != "qc" && \
      $start != "trim" && $start != "map" ]]; then
  echo "Error: problem with given argument to --start. Check the options"
  exit 1
fi

# Check that --pause is one of the possible options
if [[ $pause != "" && $pause != "indexing" && $pause != "trim" && \
      $pause != "map" && $pause != "qc" ]]; then
  echo "Error: problem with given argument to --pause. Check the options"
  exit 1
fi

# Make sure that starting step is always earlier than pausing step
if [[ $pause == "indexing" && $start != "" ]]; then
  echo "Error: problem with start, stop. Cannot stop before starting"
  exit 1
fi
if [[ $pause == "qc" && ($start != "" && $start != "qc") ]]; then
  echo "Error: problem with start, pause. Cannot stop before starting"
  exit 1
fi
if [[ $pause == "trim" && \
      ($start != "" && $start != "qc" && $start != "trim") ]]; then
  echo "Error: problem with start, pause. Cannot stop before starting"
  exit 1
fi
if [[ $start == "map" && ($pause != "" && $pause != "map") ]]; then
  echo "Error: problem with start, pause. Cannot stop before starting"
  exit 1
fi

printf "All inputs are good\n"

# Create all folders that will be needed inside output directory
if [[ $start == "" || $start == "qc" ]]; then
  rm -rf "${q}/intermediate" && \
  mkdir "${q}/intermediate" && chmod +xwr "${q}/intermediate"
  rm -rf "${q}/outputs" && \
  mkdir "${q}/outputs" && chmod +xwr "${q}/outputs"
else
  mkdir -p "${q}/intermediate" && chmod +xwr "${q}/intermediate"
  mkdir -p "${q}/outputs" && chmod +xwr "${q}/outputs"
fi

# Save number of columns in the library to use in the future
printf "num_cols:\t${ncolslib}\n" > "${q}/intermediate/useful_information.txt"

# Save input parameters that will be used in the output folder
echo "Input parameters used by CASPR" > ${q}/outputs/inputs.txt
echo "" >> ${q}/outputs/inputs.txt
echo "-q or --output-dir    = $q" >> ${q}/outputs/inputs.txt
echo "-f or --fastq-forward = $f" >> ${q}/outputs/inputs.txt
echo "-r or --fastq-reverse = $r" >> ${q}/outputs/inputs.txt
echo "-l or --library       = $l" >> ${q}/outputs/inputs.txt
echo "-e or --exper-design  = $e" >> ${q}/outputs/inputs.txt
echo "-c or --controls      = $c" >> ${q}/outputs/inputs.txt
echo "-a or --adapter-f     = $a" >> ${q}/outputs/inputs.txt
echo "-A or --adapter-r     = $A" >> ${q}/outputs/inputs.txt
echo "-m or --mismatches    = $m" >> ${q}/outputs/inputs.txt
echo "-b or --bases-aligned = $b" >> ${q}/outputs/inputs.txt
echo "-y or --fdr-threshold = $y" >> ${q}/outputs/inputs.txt
echo "-o or --orientation   = $o" >> ${q}/outputs/inputs.txt
echo "-s or --start         = $start" >> ${q}/outputs/inputs.txt
echo "-p or --pause         = $pause" >> ${q}/outputs/inputs.txt
echo "-t or --threads       = $t" >> ${q}/outputs/inputs.txt

##########
## DONE ##
##########
