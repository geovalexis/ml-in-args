# TODO: take all json files with the parsed results and create a dataframe in which:
# - each row is the sample name
# - each column is a gene and the value if whether the gene is present or not
#
# Before creating the dataframe, we need to extract all the possible genes from all the json files.


import json
import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List

import pandas as pd
import typer

logger = logging.getLogger(__name__)


@dataclass
class ResfinderResult:
    sample_name: str
    genes: List[str]
    amr_resistant: List[str]


def process_args_table(
    sample_results: ResfinderResult,
    table: Dict[str, List[int]],
    n_processed_samples: int,
):
    args_to_add = set(["_".join(arg.split("_")[:1]) for arg in sample_results.genes])
    # Add the sample values to the table
    for arg in table:
        if arg in args_to_add:
            table[arg].append(1)
            args_to_add.remove(arg)
        else:
            table[arg].append(0)
    # Add new ARG if it was not present in the previous samples
    for arg in args_to_add:
        table[arg] = [0] * (n_processed_samples - 1) + [1]


def main(
    results_files: List[str] = typer.Argument(
        ..., help="List of path to the parsed results from ResFinder"
    ),
    output: str = typer.Option(..., help="Filename of the final CSV file"),
):
    # Read and process all the results files
    samples_list = []
    args_table = defaultdict(list)
    for res_filepath in results_files:
        with open(res_filepath, "r") as json_file:
            sample_results = ResfinderResult(**json.load(json_file))
        samples_list.append(sample_results.sample_name)
        process_args_table(sample_results, args_table, len(samples_list))
    # Create dataframe
    args_table = pd.DataFrame.from_records(args_table, index=samples_list)
    args_table.to_csv(output, index_label="sample_name")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
