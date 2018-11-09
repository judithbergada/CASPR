#!/bin/bash

#####################################
## Arrange data to required format ##
#####################################

# Relationship between input parameters and the ones used here
e=$1; y=$2; q=$3; currentdir=$4; controlsfile=$5

printf "\nComputing p-values and FDRs of genes\n"

# Check the number of different controls
ctrls=$(cat $e | grep control | \
sed 's/control//g' | cut -f3 | sort -u)
# Separate different control samples according their number
i=0
for ctr in $ctrls; do
  let i=($i+1)
  c[ctr]=$(cat $e | grep "control${ctr}" | cut -f1)
done

# Check the number of different treated samples
treatm=$(cat $e | grep treated | \
sed 's/treated//g' | cut -f3 | sort -u)
# Separate different treatment samples according their number
for trm in $treatm; do
  s[trm]=$(cat $e | grep "treated${trm}" | cut -f1)
done

# Create a file with neutral controls
if [[ $controlsfile  != "" ]]; then
  cat $controlsfile | \
  awk '{if($2 == "Neutral"){print $1}}' > ${q}/intermediate/neutralctr.txt
  # Create a file with neutral guide-RNAs
  while read ctrgenename; do
    cat ${q}/outputs/table.counts.txt | \
    awk -v x="${ctrgenename}" \
    '{if($2 == x){print $1}}' >> ${q}/intermediate/neutralctrlguides.txt
  done <${q}/intermediate/neutralctr.txt
fi

############################
## Check the basic errors ##
############################

# Check if number of groups and labels are the same in control and treated
if [[ $ctrls != $treatm ]]; then
  echo "Error: problem with the experimental design file."
  echo "The numbers of control and treatment groups are different"
  exit 2
fi
# Check if there is an error with the names of the conditions
if [[ $i == 0 ]]; then
  echo "Error: problems with the experimental design file."
  echo "The format of the second column is wrong"
  exit 2
fi
# Check if SampleName's are coincident with exper.design file
allnames=$(cat $e | cut -f1)
fastqnames=$(head -n1 ${q}/outputs/table.counts.txt | cut -f 3-)
fastqnames=$(echo ",$fastqnames," | tr "[:cntrl:]" ",")
for samp in $allnames; do
  if [[ ! $fastqnames =~ ",$samp," ]]; then
    echo "Error: SampleName's of experimental design not found in fastq files"
    exit 2
  fi
done


#################################################
## Analysis of the data using MAGeCK AND PBNPA ##
#################################################

min_idx=$(echo $ctrls | awk '{print $1}')
max_idx=$(echo $ctrls | awk '{print $NF}')

for k in $ctrls; do
  ( s=$(echo ${s[$k]} | tr " +" ",")
  c=$(echo ${c[$k]} | tr " +" ",")

  # Perform the test using MAGeCK
  # If neutral controls were provided
  if [[ -s ${q}/intermediate/neutralctrlguides.txt ]]; then
    mageck test \
      -k ${q}/outputs/table.counts.txt \
      -t $s \
      -c $c \
      -n results_MAGeCK_$k \
      --control-sgrna ${q}/intermediate/neutralctrlguides.txt \
      --normcounts-to-file \
      --keep-tmp
  else # If neutral controls were not provided
    mageck test \
      -k ${q}/outputs/table.counts.txt \
      -t $s \
      -c $c \
      -n results_MAGeCK_$k \
      --normcounts-to-file \
      --keep-tmp
  fi

  # Save needed files to outputs folder
  mv results_MAGeCK_$k.gene_summary.txt \
  ${q}/outputs

  # Move the files that are not useful to intermediate folder
  mv results_MAGeCK_$k* ${q}/intermediate

  #Perfom the test using PBNPA and generate the final plots
  if [[ $k == $min_idx ]]; then
    echo "Generating plots with PBNPA results"
  fi
  Rscript --vanilla ${currentdir}/PBNPA_test.R \
          "$c" $y ${q}/outputs/table.counts.txt $e $k $q "$controlsfile" \
  || (echo "Problem with R. Check the version or the PBNPA package." && exit 2)
  if [[ $(echo $?) != 0 ]]; then exit 2; fi # Exit if there has been an error.
  if [[ $k == $max_idx ]]; then
    echo "PBNPA analysis completed successfully"
  fi ) &
done
wait

# Merge all PBNPA analysis into one pdf
gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
-sOutputFile="${q}/outputs/Selected_genes_PBNPA.pdf" \
${q}/intermediate/Selected_PBNPA*
rm ${q}/intermediate/Selected_PBNPA*

# Detect MAGeCK hits and create plots
echo "Generating plots with MAGeCK results"
Rscript --vanilla ${currentdir}/MAGeCK_test.R $y $q "$controlsfile" \
                  ${q}/intermediate/results*.gene.* \
|| (echo "Problem with R. Check the version." && exit 2)
if [[ $(echo $?) != 0 ]]; then exit 2; fi # Exit if there has been an error.
printf "MAGeCK analysis completed successfully\n"

printf "All of the analysis of genes completed succesfully\n"


##########
## DONE ##
##########
