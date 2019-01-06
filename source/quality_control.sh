#!/bin/bash

####################################
## Quality control of fastq files ##
####################################

# Relationship between input parameters and the ones used here
f=$1; r=$2; t=$3; q=$4

# Remove directory if exists and create a new one with all permissions
rm -rf "${q}/qualitycontrol" && \
mkdir "${q}/qualitycontrol" && chmod +xwr "${q}/qualitycontrol"

printf "\nPerforming the Quality Control of the fastq files\n"

# Perform quality control using FastQC tool
fastqc -t $t -o "${q}/qualitycontrol" $f
if [[ $r != "" ]]; then
  fastqc -t $t -o "${q}/qualitycontrol" $r
fi

##########
## DONE ##
##########
