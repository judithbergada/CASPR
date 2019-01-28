#!/bin/bash

set -e

#######################################
## Visualise the results using VISPR ##
#######################################

# Relationship between input parameters and the ones used here
q=$1; currentdir=$2

for res in ${q}/outputs/*.gene_summary.txt; do
  # Take the number of each test that has been performed
  testval=$(echo ${res} | \
  grep -oEi [0-9].gene.summary.txt | sed 's/.gene_summary.txt//g')
  pathname=$(echo $q | sed 's1\/1\\\/1g')

  # Create a config file for each test to display with VISPR
  cp $currentdir/config.yaml "${q}/config${testval}.yaml"

  # Modify the config file so that it contains the data for each test
  sed -i -E "s/number/${testval}/g" "${q}/config${testval}.yaml"
done

# Run VISPR even from a remote server or from a local pc
#vispr server ${q}/config*
