#!/bin/bash

# Run ResFinder4 for a certain sample
output_dir="resfinder_results"
specie="Enterococcus faecalis"
fasta_filepath="MASTER_WP4_TEAGASC_1__TGS01-01__bin.2.fa"
python -m resfinder -o $output_dir \
    -s "$specie" -l 0.6 -t 0.8 --acquired --point -ifa $fasta_filepath


