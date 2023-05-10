import logging
import os
from functools import reduce
from typing import List

import pandas as pd
import typer

logger = logging.getLogger(__name__)

column_names = ["CHROM", "POS", "REF", "ALT", "TGT", "GENE_NAME"]

NO_VARIATON_TOKEN = "N"

VARIANT_2_LABEL = {
    "A": 1,
    "G": 2,
    "C": 3,
    "T": 4,
    NO_VARIATON_TOKEN: 0,
}


def label_encode_snps(snps_table: pd.DataFrame):
    return snps_table.fillna(NO_VARIATON_TOKEN)


def process_snps_file(snps_data: pd.DataFrame, sample_name: str):
    snps_pivot = (
        snps_data[["POS", "ALT"]]
        .set_index("POS") #TODO: Set index as GENE_NAME+POSITION_IN_GENE
        .rename(columns={"ALT": sample_name})
        .transpose()
    )
    return snps_pivot


def main(
    results_files: List[str] = typer.Argument(
        ..., help="List of path to the parsed results from ResFinder"
    ),
    output: str = typer.Option(..., help="Filename of the final CSV file"),
):
    logger.info(f"Reading {len(results_files)} results files...")
    # Read and process all the results files
    snps_data_processed_all = []
    for res_file in results_files:
        sample_name = ".".join(os.path.basename(res_file).split(".")[:2])
        snps_data = pd.read_csv(res_file, sep="\t", names=column_names)
        snps_data_processed = process_snps_file(snps_data, sample_name)
        snps_data_processed.reset_index(inplace=True)
        snps_data_processed_all.append(snps_data_processed)
    # Merge all the processed SNPs
    snps_table = reduce(
        lambda left, right: pd.merge(left, right, how="outer"),
        snps_data_processed_all,
    )
    snps_table.set_index("index", inplace=True)
    snps_table = label_encode_snps(snps_table)
    # Enconde SNPs values
    snps_table = snps_table.applymap(lambda x: VARIANT_2_LABEL[x])
    # Create dataframe
    logger.info(
        f"Creating dataframe with total of {len(snps_table.columns)} different SNPs..."
    )
    snps_table.to_csv(output, index_label="sample_name")
    logger.info(f"Results saved to {output} successfully.")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
