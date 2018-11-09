#!/bin/bash

#########################
## Define basic errors ##
#########################

# Relationship between input parameters and the ones used here
e=$1; y=$2; q=$3; c=$4

printf "\nChecking that the format of inputs is correct\n"

# Check if all of the required modules are installed and
# show information on how to install the missing modules.
errorm=$(mageck test -h > /dev/null 2>&1)
if [[ $(echo $?) == 127 ]]; then
  printf "Missing: "
  echo "MAGeCK not found"
  echo "Information on the installation:"
  echo "https://bitbucket.org/liulab/mageck-vispr"
  exit 1
fi
errorm=$(vispr -h > /dev/null 2>&1)
if [[ $(echo $?) == 127 ]]; then
  printf "Missing: "
  echo "VISPR not found"
  echo "Information on the installation:"
  echo "https://bitbucket.org/liulab/mageck-vispr"
  exit 1
fi
errorm=$(R -h > /dev/null 2>&1)
if [[ $(echo $?) == 127 ]]; then
  printf "Missing: "
  echo "R not found"
  echo "Information on the installation:"
  echo "https://www.r-project.org/"
  exit 1
fi
echo "All of the required programs are properly installed"


# Check if the inputs are given in the expected formats:
# Check if the experimental design file is given as input and if it exists
if [ -z $e ]; then
  echo "Error: cannot access experimental design. Check the format of inputs"
  exit 1
fi
if [ ! -f $e ]; then
  echo "Error: the experimental design file doesn't exist"
  exit 1
fi

# Check if the controls file exists (only if it is given as input)
if [[ $c != "" ]]; then
  if [ ! -f $c ]; then
    echo "Error: the controls file doesn't exist"
    exit 1
  fi
  # Count the number of columns of the controls file
  ncolscontrol=$(head -n1 $c | awk '{print NF}')
  if [[ $ncolscontrol != 2 ]]; then
    echo "Error: the controls file doesn't have 2 columns"
    exit 1
  fi
fi

# Check if the exper.design file has hidden characters; remove them
cat $e | tr -d "\r" > ${q}/expnotused.txt
mv ${q}/expnotused.txt $e
# Count the number of columns of the experimental design
ncolsexp=$(head -n1 $e | awk '{print NF}')
if [[ $ncolsexp != 3 ]]; then
  echo "Error: the experimental design file doesn't have 3 columns"
  exit 1
fi

# Check that -y is a number between zero and 1
limits="^0\.[0-9]+$"
if [[ ! $y =~ $limits ]]; then
  echo "Error: problem with the fdr-threshold. Must be between 0 and 1"
  echo "Remember: 0 and 1 are not included."
  exit 1
fi

# Check that --pause is empty (it must be the last step if --star is test)
if [[ $pause != "" ]]; then
  echo "Error: problem with start, pause. Cannot stop before starting"
  exit 1
fi

printf "All imputs are good\n"

##########
## DONE ##
##########
