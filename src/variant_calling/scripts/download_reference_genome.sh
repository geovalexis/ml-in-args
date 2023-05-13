#!/bin/bash

# Config variables
taxon_id=$1
output_genome_path=$2
output_gff3_path=$3

# Create output directories if they don't exist
mkdir -p reference_genomes/$taxon_id

# Download reference genome and gff3 files
echo "Downloading reference genome and gff3 files for taxon $taxon_id..."
datasets download genome taxon $taxon_id --reference --assembly-source RefSeq --include genome,gff3 --filename reference_genomes/$taxon_id.zip
unzip reference_genomes/$taxon_id.zip -d reference_genomes/ncbi_$taxon_id

# Move files to outputs path
echo "Moving files to outputs path..."
current_genome_path=$(cat reference_genomes/ncbi_$taxon_id/ncbi_dataset/data/dataset_catalog.json | jq -r '.assemblies[].files[] | select(.fileType == "GENOMIC_NUCLEOTIDE_FASTA") | .filePath' | tail -n 1)
current_gff3_path=$(cat reference_genomes/ncbi_$taxon_id/ncbi_dataset/data/dataset_catalog.json | jq -r '.assemblies[].files[] | select(.fileType == "GFF3") | .filePath' | tail -n 1)
mv reference_genomes/ncbi_$taxon_id/ncbi_dataset/data/$current_genome_path $output_genome_path
mv reference_genomes/ncbi_$taxon_id/ncbi_dataset/data/$current_gff3_path $output_gff3_path

# Remove temporary files
echo "Removing temporary files..."
rm -rf reference_genomes/ncbi_$taxon_id
rm reference_genomes/$taxon_id.zip

echo "Done!"