#!/bin/bash

# Config variables
bvbrc_file=$1
output_dir=$2

# Constant variables
output_filename="PATRIC_Export.zip"

# Get list of genome IDs from BV-BRC file
genome_ids=$(awk -F "\"*,\"*" '{print $2}' $bvbrc_file | tail -n +2 | sort | uniq)

for genome_id in $genome_ids
do
  # Skip if genome already exists in output directory
  if [ -f "$output_dir/$genome_id.fna" ]; then
    echo "Genome with ID: $genome_id already exists in $output_dir"
    continue
  fi

  # Download FASTA from BV-BRC servers
  echo "Downloading Genome with ID: $genome_id..."
  curl -OJ 'https://www.bv-brc.org/api-for-website/bundle/genome/' \
    --data-raw "archiveType=zip&types=*.fna&q=in%28genome_id%2C%28$genome_id%29%29" \
    --compressed

  # Unzip and move FASTA to output directory
  unzip $output_filename
  if [ ! -d $output_dir ]; then
    mkdir $output_dir
  fi
  mv $genome_id/$genome_id.fna $output_dir

  # Clean up
  rm $output_filename
  rm -r $genome_id
done