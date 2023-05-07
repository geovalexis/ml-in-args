import json
import logging
import os
from typing import List

import pandas as pd
import typer
from lxml import etree

logger = logging.getLogger(__name__)

PHEONOTYPE_2_LABEL = {
    "intermediate": 0,
    "nonsusceptible": 0,
    "not defined": None,
    "resistant": 1,
    "susceptible": 0,
    "susceptible-dose": 0,
}


def parse_antibiogram_table(xml_text: str) -> pd.DataFrame:
    root = etree.fromstring(xml_text)
    antibiogram_table = root.xpath("//Table[@class='Antibiogram.1.0']")[0]  # type: ignore
    column_names = [col.text for col in antibiogram_table.find("Header")]  # type: ignore
    rows = [
        [cell.text for cell in row] for row in antibiogram_table.findall("Body/Row")  # type: ignore
    ]
    df = pd.DataFrame(rows, columns=column_names)
    return df


def process_amr_labels(
    antibiogram_df: pd.DataFrame,
    sample_id: str,
) -> List:
    antibiogram_df["SampleID"] = sample_id
    antibiogram_df = antibiogram_df.pivot(
        index="SampleID", columns="Antibiotic", values="Resistance phenotype"
    )
    antibiogram_df = antibiogram_df.applymap(lambda x: PHEONOTYPE_2_LABEL[x])
    antibiogram_df.reset_index(inplace=True)
    return antibiogram_df.to_dict(orient="records")


def parse_biosample_amr_labels(filepath: str) -> List:
    with open(filepath) as f:
        data = json.load(f)
    sample_uid = data["result"]["uids"][0]
    sample_accession_id = data["result"][sample_uid]["accession"]
    sample_data = data["result"][sample_uid]["sampledata"]
    logger.info(f"Parsing antibiogram table for sample '{sample_uid}'...")
    antibiogram_df = parse_antibiogram_table(sample_data)
    antibiogram_record = process_amr_labels(antibiogram_df, sample_accession_id)
    return antibiogram_record


def main(
    biosamples_summary_dir: str = typer.Option(
        ..., help="Path to the directory with biosamples summary (JSON files)"
    ),
    output: str = typer.Option(..., help="Filename of the final CSV file"),
):
    biosamples_summary_filepaths = [
        entry.path
        for entry in os.scandir(biosamples_summary_dir)
        if entry.name.endswith(".json")
    ]
    logger.info(f"Found {len(biosamples_summary_filepaths)} biosamples summary files.")
    biosamples_parsed = []
    for filepath in biosamples_summary_filepaths:
        biosamples_parsed.extend(parse_biosample_amr_labels(filepath))
    amr_labels_df = pd.DataFrame.from_records(biosamples_parsed)
    amr_labels_df.to_csv(output, index=False)
    logger.info(f"Successfully saved antibiogram table to '{output}'.")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
