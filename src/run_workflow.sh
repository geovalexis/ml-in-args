#!/bin/bash

set -e

# Config variables
config_filepath=config.yaml
cores=4

# 1. Download data (SKIP if already downloaded)
samples_dir=$(awk '/samples_dir:/ {print $2}' $config_filepath)
if [ ! -d "$samples_dir" ]; then
    echo "############## Downloading data... ##############"
    cd data_retrieval
    snakemake --use-conda --cores 4
    cd ..
fi

# 2. Perform Variant Calling analysis
cd variant_calling
echo "\n############## Performing Variant Calling... ##############"
snakemake --use-conda --cores 4
cd ..

# 3. Get ARGs genes with Resfinder
echo "\n############## Performing ARGs Calling... ##############"
cd args_calling
snakemake --use-conda --cores 4
cd ..