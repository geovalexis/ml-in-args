#!/bin/bash

# Run ResFinder4 for a certain sample
base_folder="data"
filename="1351.832.fna" # MASTER_WP4_TEAGASC_1__TGS01-01__bin.2.fa
specie="Enterococcus faecalis"
output_dir="$base_folder/resfinder_results_$filename"
fasta_filepath="$base_folder/example_inputs/$filename"
python -m resfinder -o $output_dir \
    -s "$specie" -l 0.6 -t 0.8 --acquired --point -ifa $fasta_filepath


