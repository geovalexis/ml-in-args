import logging
import os
from typing import List

import pandas as pd
import typer

logger = logging.getLogger(__name__)


def main(
    results_files: List[str] = typer.Argument(
        ..., help="List of path to the variant calling results"
    ),
    output: str = typer.Option(..., help="Filename of the final CSV file"),
):
    logger.info(f"Reading {len(results_files)} results files...")
    # Read and process all the results files
    snps_data_processed_all = []
    for res_file in results_files:
        taxid = os.path.dirname(res_file).split("/")[-1]
        sample_name = os.path.basename(res_file).rstrip(
            "filtered.snps.extract.hydrated.tsv"
        )
        snps_data = pd.read_csv(res_file, sep="\t", header=0)
        snps_data["TAX_ID"] = taxid
        snps_data["SAMPLE_NAME"] = sample_name
        snps_data_processed_all.append(snps_data)
    # Concat SNPs data
    logger.info("Concatenating SNPs data...")
    snps_table = pd.concat(snps_data_processed_all)
    snps_table.to_csv(output, index=False)
    logger.info(f"Concatenated SNPs data saved to {output}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
