#!/bin/bash

#############################################################
## Prepare merged fastq files for each of the paired reads ##
#############################################################

# Relationship between input parameters and the ones used here
f=$1; r=$2; l=$3; o=$4; a=$5; A=$6; q=$7; t=$8

printf "\nTrimming the reads\n"


############################
## Define basic functions ##
############################

# Trimming of the reads for the --orientation 35
function trimming_35 { # Only possible with forward & reverse reads
# $1=fastqfile, $2=name, $3=filename, $4=adapter-f, $5=adapter-r
# $6=lguide1, $7=lguide2, $8=ftrim, $9=fpos, $10=rtrim, $11=rpos
# Take information of first margin bases and save it to another file
  let fmargin=($9-10)
  if [[ fmargin -lt 0 ]]; then fmargin=0; fi
  let rmargin=(${11}-10)
  if [[ rmargin -lt 0 ]]; then rmargin=0; fi
  printf "${2}fmargin:\t${fmargin}\n" >> \
          "${q}/intermediate/useful_information.txt"
  printf "${2}rmargin:\t${rmargin}\n" >> \
          "${q}/intermediate/useful_information.txt"
  # Remove the first margin bases, find leftmost adapter with remaining
  # bases, trim reads and keep only as many nucleotides as length sgRNAs
  if [[ $8 == "yes" ]]; then # If trimming is needed and adapter found
    cutadapt -j $t -u ${fmargin} -g $4 -l $7 --minimum-length 5 \
    -o "${q}/intermediate/tmp.1_$2.fastq" \
    -p "${q}/intermediate/tmp.2_$2.fastq" \
    $1 $3_reverse.fastq* > "${q}/intermediate/trim_stat_$2_forward.txt"
  else # If trimming is not needed or adapter not found
    cutadapt -j $t -l $7 \
    -o "${q}/intermediate/tmp.1_$2.fastq" \
    -p "${q}/intermediate/tmp.2_$2.fastq" \
    $1 $3_reverse.fastq* > "${q}/intermediate/trim_stat_$2_forward.txt"
  fi
  echo "Forward reads finished"
  if [[ ${10} == "yes" ]]; then # Do exactly the same for reverse reads
    cutadapt -j $t -u ${rmargin} -g $5 -l $6 --minimum-length 5 \
    -o "${q}/intermediate/sgRNA1_$2.fastq" \
    -p "${q}/intermediate/sgRNA2_intermediate_$2.fastq" \
    "${q}/intermediate/tmp.2_$2.fastq" "${q}/intermediate/tmp.1_$2.fastq" > \
    "${q}/intermediate/trim_stat_$2_reverse.txt"
  else # If trimming is not needed or adapter not found
    cutadapt -j $t -l $6 \
    -o "${q}/intermediate/sgRNA1_$2.fastq" \
    -p "${q}/intermediate/sgRNA2_intermediate_$2.fastq" \
    "${q}/intermediate/tmp.2_$2.fastq" "${q}/intermediate/tmp.1_$2.fastq" > \
    "${q}/intermediate/trim_stat_$2_reverse.txt"
  fi
  echo "Reverse reads finished"
  # Remove temporal files
  rm "${q}/intermediate/tmp.1_$2.fastq" "${q}/intermediate/tmp.2_$2.fastq"
}

# Trimming of the reads for the --oreintation 53
function trimming_53 {
# $1=fastqfile, $2=name, $3=filename, $4=adapter-f, $5=adapter-r
# $6=lguide1, $7=lguide2, $8=ftrim, $9=fpos, $10=rtrim, $11=rpos
  if [[ $r != "" ]]; then # If there are forward reads
    # Take information of first margin bases and save it to another file
    let fmargin=($9-10)
    if [[ fmargin -lt 0 ]]; then fmargin=0; fi
    let rmargin=(${11}-10)
    if [[ rmargin -lt 0 ]]; then rmargin=0; fi
    printf "${2}fmargin:\t${fmargin}\n" >> \
            "${q}/intermediate/useful_information.txt"
    printf "${2}rmargin:\t${rmargin}\n" >> \
            "${q}/intermediate/useful_information.txt"
    # Remove the first margin bases, find leftmost adapter with remaining
    # bases, trim reads and keep only as many nucleotides as length sgRNAs
    if [[ $8 == "yes" ]]; then # If trimming is needed and adapter found
      cutadapt -j $t -u ${fmargin} -g $4 -l $6 --minimum-length 5 \
      -o "${q}/intermediate/tmp.1_$2.fastq" \
      -p "${q}/intermediate/tmp.2_$2.fastq" \
      $1 $3_reverse.fastq* > "${q}/intermediate/trim_stat_$2_forward.txt"
    else # If trimming is not needed or adapter not found
      cutadapt -j $t -l $6 \
      -o "${q}/intermediate/tmp.1_$2.fastq" \
      -p "${q}/intermediate/tmp.2_$2.fastq" \
      $1 $3_reverse.fastq* > "${q}/intermediate/trim_stat_$2_forward.txt"
    fi
    echo "Forward reads finished"
    if [[ ${10} == "yes" ]]; then # Do exactly the same for reverse reads
      cutadapt -j $t -u ${rmargin} -g $5 -l $7 --minimum-length 5 \
      -o "${q}/intermediate/sgRNA2_intermediate_$2.fastq" \
      -p "${q}/intermediate/sgRNA1_$2.fastq" \
      "${q}/intermediate/tmp.2_$2.fastq" "${q}/intermediate/tmp.1_$2.fastq" > \
      "${q}/intermediate/trim_stat_$2_reverse.txt"
    else # If trimming is not needed or adapter not found
      cutadapt -j $t -l $7 \
      -o "${q}/intermediate/sgRNA2_intermediate_$2.fastq" \
      -p "${q}/intermediate/sgRNA1_$2.fastq" \
      "${q}/intermediate/tmp.2_$2.fastq" "${q}/intermediate/tmp.1_$2.fastq" > \
      "${q}/intermediate/trim_stat_$2_reverse.txt"
    fi
    echo "Reverse reads finished"
    rm "${q}/intermediate/tmp.1_$2.fastq" "${q}/intermediate/tmp.2_$2.fastq"
  # If there are only forward fastq files (and not reverse)
  else
    # Take information of first margin bases and save it to another file
    let fmargin=($9-10)
    if [[ fmargin -lt 0 ]]; then fmargin=0; fi
    printf "${2}fmargin:\t${fmargin}\n" >> \
            "${q}/intermediate/useful_information.txt"
    if [[ $8 == "yes" ]]; then # If trimming is needed and adapter found
      cutadapt -j $t -u ${fmargin} -g $4 -l $6 --minimum-length 5 \
      -o "${q}/intermediate/sgRNA2_sgRNA1_$2.fastq.gz" \
      $1 > "${q}/intermediate/trim_stat_$2_forward.txt"
    else # If trimming is not needed or adapter not found
      cutadapt -j $t -l $6 \
      -o "${q}/intermediate/sgRNA2_sgRNA1_$2.fastq.gz" \
      $1 > "${q}/intermediate/trim_stat_$2_forward.txt"
    fi
    echo "Forward reads finished"
  fi
}

#####################################
## Define which trimming is needed ##
#####################################

# Automatic calculation of the gRNA lengths
lguide1=$(awk 'NR==1 {print $2}' $l | wc -c)
let lguide1=(${lguide1}-1)
if [[ $r != "" ]]; then
  lguide2=$(awk 'NR==1 {print $3}' $l | wc -c)
  let lguide2=(${lguide2}-1)
else
  lguide2=0
fi

# Automatic calculation of the adapter lengths
adaptlen=$(echo $a | wc -c)
let adaptlen=(${adaptlen}-1)
if [[ $r != "" ]]; then
  adaptlenrev=$(echo $A | wc -c)
  let adaptlenrev=(${adaptlenrev}-1)
fi

i=0
for fastqfile in $f; do
  let i=($i+1)
  echo "Analysis of ${fastqfile}:"
  # Automatic check of a repeated adapter before the guides.
  # To speed it up, consider that the first 2500 reads are representative
  # Check the most frequent position of the adapter in the reads
  echo "Checking the most frequent position of the adapter within reads"
  if file --mime-type "${fastqfile}" | grep -q gzip$; then # If compressed
    post=$(for line in $(zcat <  ${fastqfile} | awk 'NR%4==2' | \
          head -n10000 | awk 'NR%10==0'); do
            echo $line | grep -aob $a | grep -oE -m1 '[0-9]+'
          done | sort -n | uniq -c | sort -rn | head -n1)
    limitval=$(echo $post | awk '{print $2}')
    let limitval=($limitval-4)
    let upperval=($adaptlen+10)
    freq=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
    head -n10000 | awk 'NR%10==0' | \
    awk -v x="${limitval}" -v y="${upperval}" '{print substr($0, x, y)}' | \
    grep $a | wc -l)
  else # If not compressed
    post=$(for line in $(cat ${fastqfile} | awk 'NR%4==2' | \
          head -n10000 | awk 'NR%10==0'); do
            echo $line | grep -aob $a | grep -oE -m1 '[0-9]+'
          done | sort -n | uniq -c | sort -rn | head -n1)
    limitval=$(echo $post | awk '{print $2}')
    let limitval=($limitval-4)
    let upperval=($adaptlen+10)
    freq=$(cat ${fastqfile} | awk 'NR%4==2' | \
    head -n10000 | awk 'NR%10==0'| \
    awk -v x="${limitval}" -v y="${upperval}" '{print substr($0, x, y)}' | \
    grep $a | wc -l)
  fi
  val=$(echo $post | awk '{print $1}')

  # If less than 25% of adapters are found at the same position, don't trim
  if [[ $freq -lt 250 ]]; then
    trimf="no"
    positionf=0
    echo "Adapter will not be trimmed"
  else # Check that a non-constant region (the guide) follows the adapter
    echo "Checking that adapter is followed by a non-constant region (gRNA)"
    positionf=$(echo $post | awk '{print $2}')
    let positionf=($positionf + 1)
    # Check how frequent 5 nucleotides downstream the adapter are constant
    let needlen=($adaptlen + 5)
    if file --mime-type "${fastqfile}" | grep -q gzip$; then # If compressed
      cte=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
      head -n10000 | awk 'NR%10==0' | \
      awk -v x="${positionf}" -v y="${needlen}" '{print substr($0, x, y)}' | \
      grep -Ei "^$a" | awk '{print substr($0, 5, 5)}' | sort | \
      uniq -c | sort -rn | head -n1 | grep -oE "[0-9]+")
    else # If not compressed
      cte=$(cat ${fastqfile} | awk 'NR%4==2' | \
      head -n10000 | awk 'NR%10==0' | \
      awk -v x="${positionf}" -v y="${needlen}" '{print substr($0, x, y)}' | \
      grep -Ei "^$a" | awk '{print substr($0, 5, 5)}' | sort | \
      uniq -c | sort -rn | head -n1 | grep -oE "[0-9]+")
    fi

    # If the guides are found, trim. Otherwise, check second option
    if [[ $cte -lt $val/2 ]]; then
      trimf="yes"
      echo "Adapter is followed by a gRNA. It will be trimmed"
    else
      # Take the second occurence of the adapter
      let positionf=($positionf - 1)
      if file --mime-type "${fastqfile}" | grep -q gzip$; then # If compressed
        post=$(for line in \
              $(zcat <  ${fastqfile} | awk 'NR%4==2' | \
              head -n10000 | awk 'NR%10==0'); do
                echo $line | grep -aob $a | grep -oE -m2 '[0-9]+' | tail -n1
              done | \
              sed "s/${positionf}//g" | sed '/^\s*$/d' | \
              sort -n | uniq -c | sort -rn | head -n1)
        limitval=$(echo $post | awk '{print $2}')
        let limitval=($limitval-4)
        let upperval=($adaptlen+10)
        freq=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
        head -n10000 | awk 'NR%10==0' | \
        awk -v x="${limitval}" -v y="${upperval}" '{print substr($0,x,y)}' | \
        grep $a | wc -l)
      else # If not compressed
        post=$(for line in \
              $(cat ${fastqfile} | awk 'NR%4==2' | \
              head -n10000 | awk 'NR%10==0'); do
                echo $line | grep -aob $a | grep -oE -m2 '[0-9]+' | tail -n1
              done | \
              sed "s/${positionf}//g" | sed '/^\s*$/d' | \
              sort -n | uniq -c | sort -rn | head -n1)
        limitval=$(echo $post | awk '{print $2}')
        let limitval=($limitval-4)
        let upperval=($adaptlen+10)
        freq=$(cat ${fastqfile} | awk 'NR%4==2' | \
        head -n10000 | awk 'NR%10==0' | \
        awk -v x="${limitval}" -v y="${upperval}" '{print substr($0,x,y)}' | \
        grep $a | wc -l)
      fi
      val=$(echo $post | awk '{print $1}')

      # If less than 25% of adapters are found in same position, don't trim
      if [[ $freq -lt 250 ]]; then
        trimf="no"
        positionf=0
        echo "Adapter will not be trimmed. It is not followed by a gRNA"
      else # Check that a non-constant region follows the adapter
        positionf=$(echo $post | awk '{print $2}')
        let positionf=($positionf + 1)
        # Check again if nucleotides downstream the adapter are constant
        if file --mime-type "${fastqfile}" | grep -q gzip$; then # Compressed
          cte=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
          head -n10000 | awk 'NR%10==0' | awk \
          -v x="${positionf}" -v y="${needlen}" '{print substr($0, x, y)}' | \
          grep -Ei "^$a" | awk '{print substr($0, 5, 5)}' | sort | uniq -c | \
          sort -rn | head -n1 | grep -oE "[0-9]+")
        else
          cte=$(cat ${fastqfile} | awk 'NR%4==2' | \
          head -n10000 | awk 'NR%10==0' | awk \
          -v x="${positionf}" -v y="${needlen}" '{print substr($0, x, y)}' | \
          grep -Ei "^$a" | awk '{print substr($0, 5, 5)}' | sort | uniq -c | \
          sort -rn | head -n1 | grep -oE "[0-9]+")
        fi

        # If the guides are found, trim. Otherwise, don't
        if [[ $cte -lt $val/2 ]]; then
          trimf="yes"
          echo "Second occurence of adapter will be trimmed"
        else
          trimf="no"
          positionf=0
          echo "Adapter will not be trimmed. It is not followed by a gRNA"
        fi
      fi
    fi
  fi
  ftrim[$i]=$trimf
  fpos[$i]=$positionf
done

# Do exactly the same with reverse reads
if [[ ($r != "") ]]; then
  i=0
  for fastqfile in $r; do
    let i=($i+1)
    echo "Analysis of ${fastqfile}:"
    # Automatic check of a repeated adapter before the guides.
    # To speed it up, consider that the first 2500 reads are representative
    # Check the most frequent position of the adapter in the reads
    echo "Checking the most frequent position of the adapter within reads"
    if file --mime-type "${fastqfile}" | grep -q gzip$; then # If compressed
      post=$(for line in $(zcat <  ${fastqfile} | awk 'NR%4==2' | \
            head -n10000 | awk 'NR%10==0'); do
              echo $line | grep -aob $A | grep -oE -m1 '[0-9]+'
            done | sort -n | uniq -c | sort -rn | head -n1)
      limitval=$(echo $post | awk '{print $2}')
      let limitval=($limitval-4)
      let upperval=($adaptlenrev+10)
      freq=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
      head -n10000 | awk 'NR%10==0' | \
      awk -v x="${limitval}" -v y="${upperval}" '{print substr($0,x,y)}' | \
      grep $A | wc -l)
    else # If not compressed
      post=$(for line in $(cat ${fastqfile} | awk 'NR%4==2' | \
            head -n10000 | awk 'NR%10==0'); do
              echo $line | grep -aob $A | grep -oE -m1 '[0-9]+'
            done | sort -n | uniq -c | sort -rn | head -n1)
      limitval=$(echo $post | awk '{print $2}')
      let limitval=($limitval-4)
      let upperval=($adaptlenrev+10)
      freq=$(cat ${fastqfile} | awk 'NR%4==2' | \
      head -n10000 | awk 'NR%10==0' | \
      awk -v x="${limitval}" -v y="${upperval}" '{print substr($0,x,y)}' | \
      grep $A | wc -l)
    fi
    val=$(echo $post | awk '{print $1}')

    # If less than 25% of adapters are found at the same position, don't trim
    if [[ $freq -lt 250 ]]; then
      trimr="no"
      positionr=0
      echo "Adapter will not be trimmed"
    else # Check that a non-constant region (the guide) follows the adapter
      echo "Checking that adapter is followed by a non-constant region (gRNA)"
      positionr=$(echo $post | awk '{print $2}')
      let positionr=($positionr + 1)
      # Check how frequent 5 nucleotides downstream the adapter are constant
      let needlen=($adaptlenrev + 5)
      if file --mime-type "${fastqfile}" | grep -q gzip$; then # If compressed
        cte=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
        head -n10000 | awk 'NR%10==0' | \
        awk -v x="${positionr}" -v y="${needlen}" '{print substr($0,x,y)}' | \
        grep -Ei "^$A" | awk '{print substr($0,5,5)}' | \
        sort | uniq -c | sort -rn | \
        head -n1 | grep -oE "[0-9]+")
      else
        cte=$(cat ${fastqfile} | awk 'NR%4==2' | \
        head -n10000 | awk 'NR%10==0' | \
        awk -v x="${positionr}" -v y="${needlen}" '{print substr($0,x,y)}' | \
        grep -Ei "^$A" | awk '{print substr($0,5,5)}' | \
        sort | uniq -c | sort -rn | \
        head -n1 | grep -oE "[0-9]+")
      fi

      # If the guides are found, trim. Otherwise, check second option
      if [[ $cte -lt $val/2 ]]; then
        trimr="yes"
        echo "Adapter is followed by a gRNA. It will be trimmed"
      else
        # Take the second occurence of the adapter
        let positionr=($positionr - 1)
        if file --mime-type "${fastqfile}" | grep -q gzip$; then # Compressed
          post=$(for line in \
                $(zcat <  ${fastqfile} | awk 'NR%4==2' | \
                head -n10000 | awk 'NR%10==0'); do
                  echo $line | grep -aob $A | grep -oE -m2 '[0-9]+' | tail -n1
                done | \
                sed "s/${positionr}//g" | sed '/^\s*$/d' | \
                sort -n | uniq -c | sort -rn | head -n1)
          limitval=$(echo $post | awk '{print $2}')
          let limitval=($limitval-4)
          let upperval=($adaptlenrev+10)
          freq=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
          head -n10000 | awk 'NR%10==0' | \
          awk -v x="${limitval}" -v y="${upperval}" '{print substr($0,x,y)}' | \
          grep $A | wc -l)
        else # If not compressed
          post=$(for line in \
                $(cat ${fastqfile} | awk 'NR%4==2' | \
                head -n10000 | awk 'NR%10==0'); do
                  echo $line | grep -aob $A | grep -oE -m2 '[0-9]+' | tail -n1
                done | \
                sed "s/${positionr}//g" | sed '/^\s*$/d' | \
                sort -n | uniq -c | sort -rn | head -n1)
          limitval=$(echo $post | awk '{print $2}')
          let limitval=($limitval-4)
          let upperval=($adaptlenrev+10)
          freq=$(cat ${fastqfile} | awk 'NR%4==2' | \
          head -n10000 | awk 'NR%10==0' | \
          awk -v x="${limitval}" -v y="${upperval}" '{print substr($0,x,y)}' | \
          grep $A | wc -l)
        fi
        val=$(echo $post | awk '{print $1}')

        # If less than 25% of adapters are found in same position, don't trim
        if [[ $freq -lt 250 ]]; then
          trimr="no"
          positionr=0
          echo "Adapter will not be trimmed. It is not followed by a gRNA"
        else # Check that a non-constant region follows the adapter
          positionr=$(echo $post | awk '{print $2}')
          let positionr=($positionr + 1)
          # Check again if nucleotides downstream the adapter are constant
          if file --mime-type "${fastqfile}" | grep -q gzip$; then # Compressed
            cte=$(zcat <  ${fastqfile} | awk 'NR%4==2' | \
            head -n10000 | awk 'NR%10==0' | awk \
            -v x="${positionr}" -v y="${needlen}" '{print substr($0,x,y)}' | \
            grep -Ei "^$A" | awk '{print substr($0, 5, 5)}' | sort | uniq -c | \
            sort -rn | head -n1 | grep -oE "[0-9]+")
          else
            cte=$(cat ${fastqfile} | awk 'NR%4==2' | \
            head -n10000 | awk 'NR%10==0' | awk \
            -v x="${positionr}" -v y="${needlen}" '{print substr($0,x,y)}' | \
            grep -Ei "^$A" | awk '{print substr($0, 5, 5)}' | sort | uniq -c | \
            sort -rn | head -n1 | grep -oE "[0-9]+")
          fi

          # If the guides are found, trim. Otherwise, don't
          if [[ $cte -lt $val/2 ]]; then
            trimr="yes"
            echo "Second occurence of adapter will be trimmed"
          else
            trimr="no"
            positionr=0
            echo "Adapter will not be trimmed. It is not followed by a gRNA"
          fi
        fi
      fi
    fi
    rtrim[$i]=$trimr
    rpos[$i]=$positionr
  done
fi


######################
## Perform trimming ##
######################

i=0
for fastqfile in $f; do
  let i=($i+1)
  # Take only the SampleName, ignoring _reverse and directory
  name=$(echo ${fastqfile} | sed 's/_forward.fastq.*//g' | sed 's/.*\///g')
  # Remove only _forward.fastq from the name, but keep the directory info
  filename=$(echo ${fastqfile} | sed 's/_forward.fastq.*//g')

  # Do the trimming of the reads using previously defined functions
  echo "Preparing ${name} reads:"
  if [[ $o == 35 ]]; then
    trimming_35 \
    $fastqfile $name $filename $a $A $lguide1 $lguide2 \
    ${ftrim[$i]} ${fpos[$i]} ${rtrim[$i]} ${rpos[$i]}
  else
    trimming_53 \
    $fastqfile $name $filename $a $A $lguide1 $lguide2 \
    ${ftrim[$i]} ${fpos[$i]} ${rtrim[$i]} ${rpos[$i]}
  fi

  echo "Combining and preparing reads for alignment"
  # It is needed only for the paired guides.
  if [[ ($r != "") ]]; then
    # Determine encoding from fastqc outputs
    namefastqc="${q}/qualitycontrol/${name}_forward_fastqc.html"
    encoding=$(cat ${namefastqc} | grep -oEi 'Encoding.*Total')
    if [[ $(echo $encoding | grep -oEi 'Sanger') != "" || \
          $(echo $encoding | grep -oEi 'Illumina 1.9') != "" || \
          $(echo $encoding | grep -oEi 'illumina 1.8') != "" ]]; then
      phred=33
    else
      phred=64
    fi
    # Reverse-complement in order to have all reads in the same orientation
    cat "${q}/intermediate/sgRNA2_intermediate_${name}.fastq" | \
    fastx_reverse_complement -Q${phred} \
    > "${q}/intermediate/sgRNA2_${name}.fastq"
  fi

  # Paste fastq files
  if [[ $r != "" ]]; then
    paste -d"\t" \
    <(cat "${q}/intermediate/sgRNA2_${name}.fastq") \
    <(cat "${q}/intermediate/sgRNA1_${name}.fastq") | \
    awk 'BEGIN {FS="\t" }NR%2==1{print $1}NR%2==0{print $1$2}' | \
    gzip > "${q}/intermediate/sgRNA2_sgRNA1_${name}.fastq.gz"
    # Remove intermediate files that have been generated
    rm "${q}/intermediate/sgRNA1_${name}.fastq" \
        "${q}/intermediate/sgRNA2_${name}.fastq" \
        "${q}/intermediate/sgRNA2_intermediate_${name}.fastq"
  fi
  echo "${name} reads complete"
done

echo "Trimming of reads finished succesfully"

##########
## DONE ##
##########
