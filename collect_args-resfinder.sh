#!/bin/bash

# Run ResFinder4 for a certain sample
base_folder="data"
filename="1351.832.fna" # MASTER_WP4_TEAGASC_1__TGS01-01__bin.2.fa
specie="Enterococcus faecalis"
output_dir="$base_folder/resfinder_results_$filename"
fasta_filepath="$base_folder/example_inputs/$filename"
python -m resfinder -o $output_dir \
    -s "$specie" -l 0.6 -t 0.8 --acquired --point -ifa $fasta_filepath

# Results will look like as follows:
# - <fasta_filename>.json: json file with detailed results, including non-resistance genes.
# - pheno_table_species.txt: table with species specific AMR phenotypes.
# - pheno_table.txt: table with all AMR phenotypes.
# - ResFinder_Hit_in_genome_seq.fsa: fasta sequence of resistance gene hits found in the input data (query).
# - ResFinder_Resistance_gene_seq.fsa: fasta sequence of resistance gene hits found in the database (reference).
# - ResFinder_results_tab.txt: tab seperated table with predicted resistance genes.
# - ResFinder_results.txt: predicted resistance genes grouped by antibiotic class and hit alignments to reference resistance genes.
# - PointFinder_results.txt: tab seperated table with predicted point mutations leading to antibiotic resistance.
# - PointFinder_table.txt: predicted point mutations grouped into genes to which they belong.
#
# More information about the results can be found here: https://cge.food.dtu.dk/services/ResFinder-4.1/output.php
