#!/bin/bash

#################################################
## Obtain needed numbers from previous outputs ##
#################################################

# Relationship between input parameters and the ones used here
q=$1; f=$2; b=$3; l=$4; currentdir=$5

# Compute length of the guide RNAs
lguide1=$(awk 'NR==1 {print $3}' $l | wc -c)
let lguide1=(${lguide1}-1)
lguide2=$(awk 'NR==1 {print $4}' $l | wc -c)
let lguide2=(${lguide2}-1)
# Compute total length
let total_len=(${lguide1}+${lguide2}-1)


for fastqfile in $f; do
  # Take only the name of the sample file, ignoring directory
  name=$(echo ${fastqfile} | sed 's/.*\///g')

  # ________________NEEDED FOR FIRST GRAPH________________

  # Get information on number of mapped reads with user desired conditions
  mapped=$(cat $q/intermediate/Statistics_alignment_${name}.txt | \
            grep "Uniquely mapped reads number" | cut -f2)

  # Get information on total number of reads with user desired conditions
  total_reads=$(cat $q/intermediate/Statistics_alignment_${name}.txt | \
            grep "Number of input reads" | cut -f2)

  # Compute number of unmapped reads
  let unmapped=($total_reads-$mapped)

  # ________________NEEDED FOR SECOND GRAPH________________

  # Get information on number of reads that map ALL bp
  mappedm0=$(cat $q/intermediate/Statistics_alignment_${name}_totlenm0.txt | \
            grep "Uniquely mapped reads number" | cut -f2)

  # Get information on number of reads that map ALL bp with <= 3 mismatches
  mappedm3=$(cat $q/intermediate/Statistics_alignment_${name}_totlenm3.txt | \
            grep "Uniquely mapped reads number" | cut -f2)

  # Compute the number of reads that map ALL bp with more than 3 mismatches
  let mappedmore3=($mapped-$mappedm3)

  # Compute the number of reads that map ALL bp with >1 and <=3 mismatches
  let mappedm3=($mappedm3-$mappedm0)

  # Compute number of reads that are shorter than ALL bp - 3 mismatches
  let a1=($total_len+1)
  let a3=($total_len-1)
  let a4=($total_len-2)
  let a5=($total_len-3)

  short=$(cat ${q}/${name}out.sam | awk '$5 == 255' | \
          awk '{print length($10)}' | \
        grep -ve "$a1" -ve "$total_len" -ve "$a3" -ve "$a4" -ve "$a5"| wc -l)

  # Substract the number of short reads from the ones that map with > 3m
  let mappedmore3=($mappedmore3-$short)

  # ________________NEEDED FOR THIRD GRAPH________________

  # Get information on total number of reads that were unmapped previously
  reads20=$(cat $q/intermediate/Statistics_unmapped_sgrna_${name}.txt | \
            grep "Number of input reads" | cut -f2)

  # Get information on number of reads mapped to only 1 sgRNA instead of both
  mapped20=$(cat $q/intermediate/Statistics_unmapped_sgrna_${name}.txt | \
            grep "Uniquely mapped reads number" | cut -f2)

  # Get information on number of reads mapped to only 1 sgRNA multiple times
  nonunique20=$(cat $q/intermediate/Statistics_unmapped_sgrna_${name}.txt | \
            grep "Number of reads mapped to multiple loci" | cut -f2)

  # Compute number of unmapped reads
  let unmapped20=($reads20-$mapped20-$nonunique20)

  # ________________NEEDED FOR FOURTH GRAPH________________

  # Get number of multi-mapping reads with repeated guide
  repeatedg=$(cat ${q}/${name}_sgrna_out.sam | awk 'NF>5' | awk 'NF<20' | \
              awk '$12 != "NH:i:0"' | awk '$12 != "NH:i:1"' | \
              cut -f1,3 | uniq -c |  awk '$1>1' | \
              sed -r 's/(_[0-9]+)+//g' | awk '{print $2, $3}' | \
              uniq -c | wc -l)

  # Get number of multi-mapping reads mapping different genes (recombination)
  recombing=$(cat ${q}/${name}_sgrna_out.sam | awk 'NF>5' | awk 'NF<20' | \
             awk '$12 != "NH:i:0"' | awk '$12 != "NH:i:1"' | \
             cut -f1,3 | uniq -c | awk '$1==1' | \
             sed -r 's/(_[0-9]+)+//g' | awk '{print $2, $3}' | \
             uniq -c |  awk '$1==1' | awk '{print $2}' | uniq -c | wc -l)

  # Get the rest of the reads that are multi-mapped
  let others=($nonunique20-$repeatedg-$recombing)

  # ________________Plot the pie charts________________

  Rscript --vanilla ${currentdir}/alignment_pie_charts.R \
                    $name $q $b $total_len $mapped $unmapped \
                    $mappedm0 $mappedm3 $short $mappedmore3 \
                    $mapped20 $nonunique20 $unmapped20 \
                    $repeatedg $recombing $others
done

# Merge all pdf files with alignment information into one pdf
gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
-dAutoRotatePages=/None \
-sOutputFile="${q}/outputs/Alignment_statistics.pdf" \
${q}/intermediate/Alignment_stat*
rm ${q}/intermediate/Alignment_stat*

# Remove sam files
rm ${q}/*.sam

##########
## DONE ##
##########
