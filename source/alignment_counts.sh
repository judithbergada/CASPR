#!/bin/bash

#################################
## Perform useful calculations ##
#################################

# Relationship between input parameters and the ones used here
m=$1; b=$2; t=$3; q=$4; r=$5; l=$6; info=$7

# Compute length of the guide RNAs
lguide1=$(awk 'NR==1 {print $2}' $l | wc -c)
let lguide1=(${lguide1}-1)
lguide2=0
# Guide 2 will be considered only if reverse fastq are provided
if [[ $r != "" ]]; then
  lguide2=$(awk 'NR==1 {print $3}' $l | wc -c)
  let lguide2=(${lguide2}-1)
fi
let total_len=(${lguide1}+${lguide2}-1)

##############################################
## Run mapping for each of the paired reads ##
##############################################

# Remove temporal directory if exists
rm -rf "${q}/temporal"

# Perform the mapping with STAR
printf "\nAligning the reads to the guides\n"

for i in ${q}/intermediate/sgRNA2_sgRNA1*; do
  name=$(echo $i | \
  sed 's/.*intermediate\/sgRNA2_sgRNA1_//g' | sed 's/.fastq.gz//g')
  # Perform alignment
  STAR --runThreadN $t \
    --runMode alignReads \
    --genomeDir "${q}/genome" \
    --readFilesCommand zcat \
    --readFilesIn $i \
    --alignIntronMax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax $m \
    --outFilterMatchNmin $b \
    --outTmpDir "${q}/temporal" \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix "${q}/${name}_" \
    --outFilterMatchNminOverLread 0.1 \
    --outFilterMismatchNoverLmax 0.9 \
    --outFilterScoreMinOverLread 0.1 \
    --readMapNumber -1

  if [[ $info == 1 ]]; then
    # Perform alignment considering all bp
    STAR --runThreadN $t \
      --runMode alignReads \
      --genomeDir "${q}/genome" \
      --readFilesCommand zcat \
      --readFilesIn $i \
      --alignIntronMax 1 \
      --outSAMtype BAM Unsorted \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNmax 0 \
      --outFilterMatchNmin $total_len \
      --outTmpDir "${q}/temporal" \
      --outFileNamePrefix "${q}/${name}_totlenm0_" \
      --outFilterMatchNminOverLread 0.1 \
      --outFilterMismatchNoverLmax 0.9 \
      --outFilterScoreMinOverLread 0.1 \
      --readMapNumber -1

      # Perform alignment considering all bp and 3 mismatches
      STAR --runThreadN $t \
        --runMode alignReads \
        --genomeDir "${q}/genome" \
        --readFilesCommand zcat \
        --readFilesIn $i \
        --alignIntronMax 1 \
        --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 3 \
        --outFilterMatchNmin $total_len \
        --outTmpDir "${q}/temporal" \
        --outFileNamePrefix "${q}/${name}_totlenm3_" \
        --outFilterMatchNminOverLread 0.1 \
        --outFilterMismatchNoverLmax 0.9 \
        --outFilterScoreMinOverLread 0.1 \
        --readMapNumber -1

      # Save useful files considering all bp
      samtools view -h -o "${q}/${name}out.sam" \
                          "${q}/${name}_Aligned.sortedByCoord.out.bam"
      mv "${q}/${name}_totlenm0_Log.final.out" \
      ${q}/intermediate/Statistics_alignment_${name}_totlenm0.txt
      mv "${q}/${name}_totlenm3_Log.final.out" \
      ${q}/intermediate/Statistics_alignment_${name}_totlenm3.txt

      # Remove files that can lead to problems in future iterations
      rm ${q}/${name}_totlenm0_Aligned* \
          ${q}/${name}_totlenm0_Log* ${q}/${name}_totlenm0_SJ.out.tab
      rm ${q}/${name}_totlenm3_Aligned* \
          ${q}/${name}_totlenm3_Log* ${q}/${name}_totlenm3_SJ.out.tab
  fi

  # Compute the number of reads per sample
  samtools index "${q}/${name}_Aligned.sortedByCoord.out.bam"
  samtools idxstats "${q}/${name}_Aligned.sortedByCoord.out.bam" | \
  cut -f1,3 > "${q}/intermediate/reads_${name}.tsv"

  # Save useful files considering user parameters
  mv "${q}/${name}_Log.final.out" \
  ${q}/intermediate/Statistics_alignment_${name}.txt
  mv "${q}/${name}_Unmapped.out.mate1" \
  ${q}/intermediate/Unmapped_${name}.fastq

  # Remove files that can lead to problems in future iterations
  rm ${q}/${name}_Aligned* ${q}/${name}_Log* ${q}/${name}_SJ.out.tab

done


###################################################
## Prepare final table with the counts per guide ##
###################################################

# Write first column of the table: IDs of the guide RNAs
commandpaste="<(sort -V ${q}/intermediate/reads_${name}.tsv | cut -f1)"

# Write second column of the table: Gene name
commandpaste="${commandpaste} \
<(sort -V ${q}/intermediate/reads_${name}.tsv | cut -f1 | sed 's/_.*//g')"

# Write header
header="ID\tGene"

# Write the rest of the columns: number of counts per SampleName
for i in ${q}/intermediate/reads*; do
  commandpaste="${commandpaste} <(sort -V $i | cut -f2)"
  name=$(echo $i | sed 's/.*intermediate\/reads_//g' | sed 's/.tsv//g')
  header="${header}\t${name}"
done

# Paste everything
eval paste ${commandpaste} > "${q}/intermediate/counts.txt"
echo -e ${header} | \
cat - ${q}/intermediate/counts.txt > "${q}/outputs/table.counts.txt"


###################################################################
## Check what happens with unmapped reads only for paired-guides ##
###################################################################

if [[ $info == 1 ]]; then

  # If the pgRNAs have different length, take the lenght of the shortest
  lenguide=$( [ ${lguide1} -lt ${lguide2} ] \
        && echo "${lguide1}" || echo "${lguide2}" )

  # Perform the mapping of the unmapped reads with STAR
  printf "\nChecking unmapped reads\n"

  # Check if at least one of the guide RNAs (the shortest) is aligned
  for i in ${q}/intermediate/Unmapped_*; do
    name=$(echo $i | \
    sed 's/.*intermediate\/Unmapped_//g' | sed 's/.fastq//g')
    # Perform alignment
    STAR --runThreadN $t \
      --runMode alignReads \
      --genomeDir "${q}/genome" \
      --readFilesIn $i \
      --alignIntronMax 1 \
      --outSAMunmapped Within \
      --outSAMtype BAM Unsorted \
      --outFilterMultimapNmax 20 \
      --outReadsUnmapped Fastx \
      --outFilterMismatchNmax 0 \
      --outFilterMatchNmin $lenguide \
      --outTmpDir "${q}/temporal" \
      --outFileNamePrefix "${q}/${name}_sgrna_" \
      --outFilterMatchNminOverLread 0.1 \
      --outFilterMismatchNoverLmax 0.9 \
      --outFilterScoreMinOverLread 0.1 \
      --readMapNumber -1

    # Save statistics
    mv "${q}/${name}_sgrna_Log.final.out" \
    ${q}/intermediate/Statistics_unmapped_sgrna_${name}.txt
    mv "${q}/${name}_sgrna_Unmapped.out.mate1" \
    ${q}/intermediate/Unmapped_sgrna_${name}.fastq
    samtools view -h -o "${q}/${name}_sgrna_out.sam" \
                        "${q}/${name}_sgrna_Aligned.out.bam"
    # Remove files that can lead to problems in future iterations
    rm ${q}/${name}_sgrna_Aligned* \
        ${q}/${name}_sgrna_Log* ${q}/${name}_sgrna_SJ.out.tab
  done
  rm ${q}/intermediate/Unmapped_*
fi

# Remove any file introduced by aligner and not needed
rm -f ${q}/temporal*
rm -rf ${q}/temporal*

##########
## DONE ##
##########
