import logging

import numpy as np
import pandas as pd
import typer

logger = logging.getLogger(__name__)

MEASUREMENT_SIGN_2_PHENOTYPE = {
    "<": "Susceptible",
    "<=": "Susceptible",
    ">=": "Resistant",
    "=": "Susceptible",
    ">": "Resistant",
    "==": "Intermediate",
    np.nan: np.nan,
}


def process_amr_labels(
    bvbrc_df: pd.DataFrame,
) -> pd.DataFrame:
    # TODO: move this logic to data cleaning step (in R)
    bvbrc_df_no_dups = bvbrc_df.drop_duplicates()
    bvbrc_df_no_dups = bvbrc_df_no_dups.drop_duplicates(["Genome ID", "Antibiotic"])
    bvbrc_df_no_dups.where(bvbrc_df_no_dups["Resistant Phenotype"] != np.nan)
    bvbrc_df_no_dups["Resistant Phenotype"] = bvbrc_df_no_dups.apply(
        lambda x: MEASUREMENT_SIGN_2_PHENOTYPE[x["Measurement Sign"]]
        if pd.isna(x["Resistant Phenotype"])
        else x["Resistant Phenotype"],
        axis=1,
    )
    bvbrc_pivot = bvbrc_df_no_dups.pivot(
        index="Genome ID", columns="Antibiotic", values="Resistant Phenotype"
    )
    return bvbrc_pivot


def main(
    bvbrc_file: str = typer.Option(..., help="Path to the BVBRC CSV file"),
    output: str = typer.Option(..., help="Filename of the final CSV file"),
):
    logger.info(f"Reading BVBRC file from '{bvbrc_file}'...")
    # Read and process file
    bvbrc_df = pd.read_csv(bvbrc_file, header=0)
    bvbrc_df_processed = process_amr_labels(bvbrc_df)
    # Save file
    logger.info(f"Processed file has the following shape: {bvbrc_df_processed.shape}")
    bvbrc_df_processed.to_csv(output, index=True, header=True)
    logger.info(f"Processed file successfully saved to {output}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
