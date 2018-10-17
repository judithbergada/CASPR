#!/bin/bash

##############################################
## Run mapping for each of the paired reads ##
##############################################

# Relationship between input parameters and the ones used here
m=$1; b=$2; t=$3; q=$4; r=$5; l=$6

# Remove temporal directory if exists
rm -rf "${q}temporal"

# Perform the mapping with STAR
printf "\nAligning the reads to the guides\n"

for i in ${q}intermediate/sgRNA2_sgRNA1*; do
  name=$(echo $i | \
  sed 's/.*intermediate\/sgRNA2_sgRNA1_//g' | sed 's/.fastq.gz//g')
  # Perform alignment
  STAR --runThreadN $t \
  --runMode alignReads \
  --genomeDir "${q}genome" \
  --readFilesCommand zcat \
  --readFilesIn $i \
  --alignIntronMax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNmax $m \
  --outFilterMatchNmin $b \
  --outTmpDir "${q}temporal" \
  --outReadsUnmapped Fastx \
  --outFileNamePrefix "${q}${name}_" \
  --outFilterMatchNminOverLread 0.1 \
  --outFilterMismatchNoverLmax 0.9 \
  --outFilterScoreMinOverLread 0.1 \
  --readMapNumber -1

  # Compute the number of reads per sample
  samtools index "${q}${name}_Aligned.sortedByCoord.out.bam"
  samtools idxstats "${q}${name}_Aligned.sortedByCoord.out.bam" | \
  cut -f1,3 > "${q}intermediate/reads_${name}.tsv"

  # Save useful files
  mv "${q}${name}_Log.final.out" \
  ${q}intermediate/Statistics_alignment_${name}.txt
  mv "${q}${name}_Unmapped.out.mate1" \
  ${q}intermediate/Unmapped_${name}.fastq

  # Remove files that can lead to problems in future iterations
  rm ${q}${name}_Aligned* ${q}${name}_Log* ${q}${name}_SJ.out.tab
done


###################################################
## Prepare final table with the counts per guide ##
###################################################

# Write first column of the table: IDs of the guide RNAs
commandpaste="<(sort -V ${q}intermediate/reads_${name}.tsv | cut -f1)"

# Write second column of the table: Gene name
commandpaste="${commandpaste} \
<(sort -V ${q}intermediate/reads_${name}.tsv | cut -f1 | sed 's/_.*//g')"

# Write header
header="ID\tGene"

# Write the rest of the columns: number of counts per SampleName
for i in ${q}intermediate/reads*; do
  commandpaste="${commandpaste} <(sort -V $i | cut -f2)"
  name=$(echo $i | sed 's/.*intermediate\/reads_//g' | sed 's/.tsv//g')
  header="${header}\t${name}"
done

# Paste everything
eval paste ${commandpaste} > "${q}intermediate/counts.txt"
echo -e ${header} | \
cat - ${q}intermediate/counts.txt > "${q}outputs/table.counts.txt"


###################################################################
## Check what happens with unmapped reads only for paired-guides ##
###################################################################

if [[ $r != "" ]]; then

  # Compute length of the guide RNAs
  lguide1=$(awk 'NR==1 {print $2}' $l | wc -c)
  lguide2=$(awk 'NR==1 {print $3}' $l | wc -c)
  # If the pgRNAs have different length, take the lenght of the shortest
  lenguide=$( [ ${lguide1} -lt ${lguide2} ] \
        && echo "${lguide1}" || echo "${lguide2}" )
  let lenguide=(${lenguide}-1)
  # Compute an alternative length longer than lenguide
  let lenaltern=(${lenguide}+5)

  # Perform the mapping of the unmapped reads with STAR
  printf "\nChecking unmapped reads\n"

  # Check if at least one of the guide RNAs (the shortest) is aligned
  for i in ${q}intermediate/Unmapped_*; do
    name=$(echo $i | \
    sed 's/.*intermediate\/Unmapped_//g' | sed 's/.fastq//g')
    # Perform alignment
    STAR --runThreadN $t \
    --runMode alignReads \
    --genomeDir "${q}genome" \
    --readFilesIn $i \
    --alignIntronMax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax $m \
    --outFilterMatchNmin $lenguide \
    --outTmpDir "${q}temporal" \
    --outFileNamePrefix "${q}${name}_" \
    --outFilterMatchNminOverLread 0.1 \
    --outFilterMismatchNoverLmax 0.9 \
    --outFilterScoreMinOverLread 0.1 \
    --readMapNumber -1

    # Save statistics
    mv "${q}${name}_Log.final.out" \
    ${q}intermediate/Statistics_unmapped_sgrna_${name}.txt
    # Remove files that can lead to problems in future iterations
    rm ${q}${name}_Aligned* ${q}${name}_Log* ${q}${name}_SJ.out.tab

    # Check if part of the second cut is also aligned.
    # In this case, one of the gRNAs has been cut during the sequencing.
    STAR --runThreadN $t \
    --runMode alignReads \
    --genomeDir "${q}genome" \
    --readFilesIn $i \
    --alignIntronMax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax $m \
    --outFilterMatchNmin $lenaltern \
    --outTmpDir "${q}temporal" \
    --outFileNamePrefix "${q}${name}_" \
    --outFilterMatchNminOverLread 0.1 \
    --outFilterMismatchNoverLmax 0.9 \
    --outFilterScoreMinOverLread 0.1 \
    --readMapNumber -1

      # Save statistics
      mv "${q}${name}_Log.final.out" \
      ${q}intermediate/Statistics_unmapped_pgrna_${name}.txt
      # Remove files that can lead to problems in future iterations
      rm ${q}${name}_Aligned* ${q}${name}_Log* ${q}${name}_SJ.out.tab
  done
fi

# Remove fastq files with unmapped reads
#rm ${q}intermediate/Unmapped*

##########
## DONE ##
##########
