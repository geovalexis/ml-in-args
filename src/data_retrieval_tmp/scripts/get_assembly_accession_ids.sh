#!/bin/bash

# Config variables
input_filepath=$1
output_filepath=$2


# Read in list of biosample accession IDs
biosamples_accession_ids=$(awk -F , '{print $1}' "$input_filepath" | tail -n +2 | tail -n 2)

# Load and prepare progress bar
source ./scripts/progress_bar.sh
i=1

# Get assembly accession IDs from biosample accession IDs
n_samples=$(echo $biosamples_accession_ids | wc -w)
for id in $biosamples_accession_ids
do
    esearch -db biosample -query $id | elink -target assembly | esummary | xtract \
    -pattern DocumentSummary -element AssemblyAccession >> "$output_filepath"
    ProgressBar $i $n_samples
    i=$((i+1))
done