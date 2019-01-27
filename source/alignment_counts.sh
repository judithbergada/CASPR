#!/bin/bash

set -e

#################################
## Perform useful calculations ##
#################################

# Relationship between input parameters and the ones used here
m=$1; b=$2; t=$3; q=$4; r=$5; l=$6; info=$7

# Compute length of the guide RNAs
lguide1=$(awk 'NR==1 {print $3}' $l | wc -c)
let lguide1="${lguide1}-1"
lguide2=0
# Guide 2 will be considered only if reverse fastq are provided
if [[ $r != "" ]]; then
  lguide2=$(awk 'NR==1 {print $4}' $l | wc -c)
  let lguide2="${lguide2}-1"
fi
let total_len="${lguide1}+${lguide2}-1"

# Check the Shared Memory available
available_shm=$(sysctl -A 2>/dev/null | grep shmmax | grep -Eo "[0-9]+")
req_shm=10000000000 # 10GB
if [ "$available_shm" -lt "$req_shm" ]; then
  shm_flags=""
else
  shm_flags="--genomeLoad LoadAndKeep --limitBAMsortRAM $req_shm"
  echo "Using shared memory to load the genome."
fi

##############################################
## Run mapping for each of the paired reads ##
##############################################

# Remove temporal directory if exists
rm -rf "${q}/temporal"

# Perform the mapping with STAR
printf "\nAligning the reads to the guides\n"

for i in ${q}/intermediate/sgRNA2_sgRNA1*; do
  # Take name of fastqfile ignoring directory, zip, fastq and common part
  nametwo=$(echo $i | sed 's/.*intermediate\/sgRNA2_sgRNA1_//g' | \
        sed 's/\.gz//g' | sed 's/\.fastq//g' | sed 's/\.fq//g')

  # Perform alignment
  STAR --runThreadN $t \
    --runMode alignReads \
    --genomeDir "${q}/genome" \
    --readFilesCommand "gunzip -c" \
    --readFilesIn $i \
    $shm_flags \
    --alignIntronMax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax $m \
    --outFilterMatchNmin $b \
    --outTmpDir "${q}/temporal" \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix "${q}/${nametwo}_" \
    --outFilterMatchNminOverLread 0.1 \
    --outFilterMismatchNoverLmax 0.9 \
    --outFilterScoreMinOverLread 0.1 \
    --readMapNumber -1

  if [[ $info == 1 ]]; then
    # Perform alignment considering all bp
    STAR --runThreadN $t \
      --runMode alignReads \
      --genomeDir "${q}/genome" \
      --readFilesCommand "gunzip -c" \
      --readFilesIn $i \
      $shm_flags \
      --alignIntronMax 1 \
      --outSAMtype BAM Unsorted \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNmax 0 \
      --outFilterMatchNmin $total_len \
      --outTmpDir "${q}/temporal" \
      --outFileNamePrefix "${q}/${nametwo}_totlenm0_" \
      --outFilterMatchNminOverLread 0.1 \
      --outFilterMismatchNoverLmax 0.9 \
      --outFilterScoreMinOverLread 0.1 \
      --readMapNumber -1

      # Perform alignment considering all bp and 3 mismatches
      STAR --runThreadN $t \
        --runMode alignReads \
        --genomeDir "${q}/genome" \
        --readFilesCommand "gunzip -c" \
        --readFilesIn $i \
        $shm_flags \
        --alignIntronMax 1 \
        --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 3 \
        --outFilterMatchNmin $total_len \
        --outTmpDir "${q}/temporal" \
        --outFileNamePrefix "${q}/${nametwo}_totlenm3_" \
        --outFilterMatchNminOverLread 0.1 \
        --outFilterMismatchNoverLmax 0.9 \
        --outFilterScoreMinOverLread 0.1 \
        --readMapNumber -1

      # Save useful files considering all bp
      samtools view -h -o "${q}/${nametwo}out.sam" \
                          "${q}/${nametwo}_Aligned.sortedByCoord.out.bam"
      mv "${q}/${nametwo}_totlenm0_Log.final.out" \
      ${q}/intermediate/Statistics_alignment_${nametwo}_totlenm0.txt
      mv "${q}/${nametwo}_totlenm3_Log.final.out" \
      ${q}/intermediate/Statistics_alignment_${nametwo}_totlenm3.txt

      # Remove files that can lead to problems in future iterations
      rm ${q}/${nametwo}_totlenm0_Aligned* \
          ${q}/${nametwo}_totlenm0_Log* ${q}/${nametwo}_totlenm0_SJ.out.tab
      rm ${q}/${nametwo}_totlenm3_Aligned* \
          ${q}/${nametwo}_totlenm3_Log* ${q}/${nametwo}_totlenm3_SJ.out.tab
  fi

  # Compute the number of reads per sample
  samtools index "${q}/${nametwo}_Aligned.sortedByCoord.out.bam"
  samtools idxstats "${q}/${nametwo}_Aligned.sortedByCoord.out.bam" | \
  cut -f1,3 > "${q}/intermediate/reads_${nametwo}.tsv"

  # Remove rows that have been added at the end without information
  cat "${q}/intermediate/reads_${nametwo}.tsv" | \
        grep -v -E "^\*\s0$" > "${q}/intermediate/reads_new.tsv"
  mv "${q}/intermediate/reads_new.tsv" "${q}/intermediate/reads_${nametwo}.tsv"

  # Save useful files considering user parameters
  mv "${q}/${nametwo}_Log.final.out" \
  ${q}/intermediate/Statistics_alignment_${nametwo}.txt
  mv "${q}/${nametwo}_Unmapped.out.mate1" \
  ${q}/intermediate/Unmapped_${nametwo}

  # Remove files that can lead to problems in future iterations
  rm ${q}/${nametwo}_Aligned* ${q}/${nametwo}_Log* ${q}/${nametwo}_SJ.out.tab

done


###################################################
## Prepare final table with the counts per guide ##
###################################################

# Write first column of the table: IDs of the guide RNAs
commandpaste="<(sort -V ${q}/intermediate/reads_${nametwo}.tsv | cut -f1)"

# Write second column of the table: Gene name
commandpaste="${commandpaste} \
    <(sort -u -k3 ${q}/intermediate/sgRNA2.sgRNA1_map.txt | sort -V | cut -f2)"

# Write header
header="ID\tGene"

# Write the rest of the columns: number of counts per SampleName
for i in ${q}/intermediate/reads*; do
  commandpaste="${commandpaste} <(sort -V $i | cut -f2)"
  name=$(echo $i | sed 's/.*intermediate\/reads_//g' | sed 's/.tsv//g')
  nametwo=$(echo ${q}/intermediate/sgRNA2_sgRNA1_${name}* | \
          sed 's/.*intermediate\/sgRNA2_sgRNA1_//g')
  header="${header}\t${nametwo}"
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
    # Take name ignoring directory and format
    nametwo=$(echo $i | sed 's/.*intermediate\/Unmapped_//g' | \
          sed 's/\.gz//g' | sed 's/\.fastq//g' | sed 's/\.fq//g')

    # Perform alignment
    STAR --runThreadN $t \
      --runMode alignReads \
      --genomeDir "${q}/genome" \
      --readFilesIn $i \
      $shm_flags \
      --alignIntronMax 1 \
      --outSAMunmapped Within \
      --outSAMtype BAM Unsorted \
      --outFilterMultimapNmax 20 \
      --outReadsUnmapped Fastx \
      --outFilterMismatchNmax 0 \
      --outFilterMatchNmin $lenguide \
      --outTmpDir "${q}/temporal" \
      --outFileNamePrefix "${q}/${nametwo}_sgrna_" \
      --outFilterMatchNminOverLread 0.1 \
      --outFilterMismatchNoverLmax 0.9 \
      --outFilterScoreMinOverLread 0.1 \
      --readMapNumber -1

    # Save statistics
    mv "${q}/${nametwo}_sgrna_Log.final.out" \
    ${q}/intermediate/Statistics_unmapped_sgrna_${nametwo}.txt
    mv "${q}/${nametwo}_sgrna_Unmapped.out.mate1" \
    ${q}/intermediate/Unmapped_sgrna_${nametwo}
    samtools view -h -o "${q}/${nametwo}_sgrna_out.sam" \
                        "${q}/${nametwo}_sgrna_Aligned.out.bam"
    # Remove files that can lead to problems in future iterations
    rm ${q}/${nametwo}_sgrna_Aligned* \
        ${q}/${nametwo}_sgrna_Log* ${q}/${nametwo}_sgrna_SJ.out.tab
  done
  rm ${q}/intermediate/Unmapped_*
fi

# Remove any file introduced by aligner and not needed
if [[ $shm_flags != "" ]]; then
  STAR --runThreadN $t \
    --runMode alignReads \
    --genomeDir "${q}/genome" \
    --genomeLoad Remove \
    --outTmpDir "${q}/temporal" \
    --outFileNamePrefix "${q}/removalprocess" \
    --limitBAMsortRAM $req_shm
    rm ${q}/removalprocess*
fi
rm -rf ${q}/temporal*
rm -f ${q}/temporal*

##########
## DONE ##
##########
