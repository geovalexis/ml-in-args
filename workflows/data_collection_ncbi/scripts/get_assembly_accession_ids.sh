#!/bin/bash

# Config variables
input_filepath=$1
output_filepath=$2

# Read in list of biosample accession IDs
biosamples_accession_ids=$(awk -F , '{print $1}' "$input_filepath" | tail -n +2)

# Create query to get assembly accession IDs from biosample accession IDs
echo "Getting assembly accessionIds for total of $(echo $biosamples_accession_ids | wc -w) biosamples..."
queryIds=$(echo $biosamples_accession_ids | sed 's/ / OR /g')
esearch -db biosample -query "$queryIds" | elink -target assembly | esummary | xtract \
-pattern DocumentSummary -element BioSampleAccn,AssemblyAccession,SpeciesTaxid > "$output_filepath"
echo "Done. Output written to $output_filepath"