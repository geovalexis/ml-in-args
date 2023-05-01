import logging

import pandas as pd
import typer

logger = logging.getLogger(__name__)

SIGN_2_LABEL = {
    "<=": 0,
    "=": 0,
    ">": 1,
}


def process_amr_labels(
    bvbrc_df: pd.DataFrame,
) -> pd.DataFrame:
    bvbrc_pivot = bvbrc_df.pivot(
        index="Genome ID", columns="Antibiotic", values="Measurement Sign"
    )
    bvbrc_pivot = bvbrc_pivot.applymap(lambda x: SIGN_2_LABEL[x])
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
