#!/bin/bash

########################
## Create fake genome ##
########################

# Relationship between input parameters and the ones used here
l=$1; t=$2; q=$3

# Remove directory if exists and create a new one with all permissions
rm -rf "${q}genome" && mkdir "${q}genome" && chmod +xwr "${q}genome"

printf "\nCreating fake genome with the guide RNAs\n"

# Create fake genome.fasta file
# Check if sgRNAs or pgRNAs
ncolslib=$(head -n1 "${q}intermediate/useful_information.txt" | \
            awk '{print $2}')
if [[ $ncolslib == 3 ]]; then
  cat $l | awk '{print $1, $3 $2}' > "${q}intermediate/sgRNA2.sgRNA1_map.txt"
else
  cp $l "${q}intermediate/sgRNA2.sgRNA1_map.txt"
fi
# Transform the genome to a fasta format
sort -u -k2 "${q}intermediate/sgRNA2.sgRNA1_map.txt" | \
sed 's/^[\t]*/>/g;s/\s/\n/g' > "${q}genome/genome.fasta"

# Index genome
STAR --runThreadN $t \
--runMode genomeGenerate \
--genomeSAindexNbases 8 \
--outFileNamePrefix "${q}genome/" \
--outTmpDir "${q}temporal" \
--genomeDir "${q}genome" \
--genomeFastaFiles "${q}genome/genome.fasta"


##########
## DONE ##
##########
