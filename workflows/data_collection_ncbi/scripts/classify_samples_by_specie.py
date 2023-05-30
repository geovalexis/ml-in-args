import logging
import os

import pandas as pd
import typer

logger = logging.getLogger(__name__)


def main(
    accession_2_taxid: str = typer.Option(
        ..., help="Tabular file with accession to taxid mapping"
    ),
    samples_dir: str = typer.Option(..., help="Name of the directory with the samples"),
):
    samples = [
        sample for sample in os.scandir(samples_dir) if sample.name.endswith(".fna")
    ]
    accession_2_taxid_df = pd.read_csv(
        accession_2_taxid, sep="\t", names=["accession", "taxid"]
    )
    logger.info("Classifying total of %d samples...", len(samples))
    for sample in samples:
        sample_accession = "_".join(sample.name.split("_")[0:2])
        try:
            sample_taxid = accession_2_taxid_df[
                sample_accession == accession_2_taxid_df.accession
            ].taxid.values[0]
            sample_dir = os.path.join(samples_dir, str(sample_taxid))
            if not os.path.isdir(sample_dir):
                os.mkdir(sample_dir)
            os.rename(sample.path, os.path.join(sample_dir, sample.name))
        except IndexError:
            logger.warning(f"Sample '{sample_accession}' not found in accession file")
            continue
    logger.info("Done!")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
