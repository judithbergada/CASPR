#!/bin/bash

########################
## Create fake genome ##
########################

# Relationship between input parameters and the ones used here
l=$1; t=$2; q=$3

# Remove directory if exists and create a new one with all permissions
rm -rf "${q}/genome" && mkdir "${q}/genome" && chmod +xwr "${q}/genome"

printf "\nCreating fake genome with the guide RNAs\n"

# Create fake genome.fasta file
# Check if sgRNAs or pgRNAs
ncolslib=$(head -n1 "${q}/intermediate/useful_information.txt" | \
            awk '{print $2}')
if [[ $ncolslib == 4 ]]; then
  cat $l | awk -v OFS='\t' '{print $1, $2, $4 $3}' \
            > "${q}/intermediate/sgRNA2.sgRNA1_map.txt"
else
  cat $l | awk -v OFS='\t' '{print $1, $2, $3}' \
            > "${q}/intermediate/sgRNA2.sgRNA1_map.txt"
fi
# Transform the genome to a fasta format
sort -u -k3 "${q}/intermediate/sgRNA2.sgRNA1_map.txt" | \
awk -v OFS='\t' '{print $1, $3}' | \
sed 's/^[\t]*/>/g;s/\s/\n/g' > "${q}/genome/genome.fasta"

# Compute needed --genomeSaindexNbases parameter according to STAR formula
# Formula: genomeSaindexNbases = log2(numbases)/2 - 1
num_bases_grna=$(cat ${q}/genome/genome.fasta | awk 'NR==4' | wc -c)
num_grnas=$(wc -l ${q}/genome/genome.fasta | awk '{print $1/2}')
let totalbp=(${num_bases_grna}*${num_grnas})
log2res=$(echo "l($totalbp)/l(2)" | bc -l)
genomeSaind=$(echo ${log2res} | awk '{printf "%.0f\n", $1/2 - 2}')

# Index genome
STAR --runThreadN $t \
--runMode genomeGenerate \
--genomeSAindexNbases ${genomeSaind} \
--outFileNamePrefix "${q}/genome/" \
--outTmpDir "${q}/temporal" \
--genomeDir "${q}/genome" \
--genomeFastaFiles "${q}/genome/genome.fasta"

##########
## DONE ##
##########
