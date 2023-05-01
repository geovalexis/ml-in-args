import json
import logging
from dataclasses import dataclass
from typing import Any, Dict, List

import typer

logger = logging.getLogger(__name__)

# Resfinder results will be saved into a directory with the sample's name and it will look like as follows:
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


@dataclass
class SeqRegion:
    type: str
    phenotypes: List[str]
    ref_database: List[str]
    gene: bool
    ref_id: str
    name: str
    ref_acc: str
    identity: float
    alignment_length: int
    ref_seq_lenght: int
    ref_start_pos: int
    ref_end_pos: int
    query_id: str
    query_start_pos: int
    query_end_pos: int
    pmids: List[int]
    notes: List[str]
    coverage: int
    grade: int
    key: str


@dataclass
class Phenotype:
    type: str
    amr_classes: List[str]
    seq_regions: List[str]
    seq_variations: List[Any]
    ref_database: List[str]
    category: str
    key: str
    amr_resistance: str
    amr_resistant: bool
    amr_species_relevant: bool
    grade: int


def main(
    results_dir: str = typer.Option(..., help="Path to the results folder from ResFinder"),
    output: str = typer.Option(..., help="Path to the output file"),
):
    sample_name = results_dir.split("/")[-1]
    # Read json file with results
    logger.info(f"Reading results for sample {sample_name}...")
    json_results_filepath = f"{results_dir}/{sample_name.split('.')[0]}.json"
    with open(json_results_filepath, "r") as json_file:
        resfinder_results: Dict = json.load(json_file)
    # Parse only those sections that are relevant for our analysis
    genes = [SeqRegion(**seq_region) for _, seq_region in resfinder_results["seq_regions"].items()]
    logger.info(f"Found {len(genes)} genes.")
    phenotypes = [Phenotype(**phenotype) for _, phenotype in resfinder_results["phenotypes"].items()]
    resistant_antiobiotics = list(filter(lambda x: x.amr_resistant, phenotypes))
    logger.info(f"Found {len(resistant_antiobiotics)}/{len(phenotypes)} antibiotics with resistance")
    # Create output file
    output_data = {
        "sample_name": sample_name,
        "genes": [gene.ref_id for gene in genes],
        "amr_resistant": [antibiotic.amr_resistance for antibiotic in resistant_antiobiotics],
    }
    with open(output, "w") as output_file:
        json.dump(output_data, output_file, indent=4)
    logger.info(f"Results saved to {output} successfully.")

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
