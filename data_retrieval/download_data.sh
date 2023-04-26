#!/bin/bash

# Config variables
bvbrc_file=$1
output_filename="PATRIC_Export.zip"

genome_ids=$(awk -F "\"*,\"*" '{print $2}' $bvbrc_file | tail -n +2)
for genome_id in $genome_ids
do
  echo "Downloading Genome with ID: $genome_id..."

  # Download FASTA from BV-BRC servers
  curl -OJ 'https://www.bv-brc.org/api-for-website/bundle/genome/' \
    --data-raw "archiveType=zip&types=*.fna&q=in%28genome_id%2C%28$genome_id%29%29" \
    --compressed

  # Unzip and move FASTA to data/samples
  unzip $output_filename
  mv $genome_id/$genome_id.fna data/samples

  # Clean up
  rm $output_filename
  rm -r $genome_id
done