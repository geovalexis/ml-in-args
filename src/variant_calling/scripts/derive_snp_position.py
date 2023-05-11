import logging

import pandas as pd
import typer

logger = logging.getLogger(__name__)

REF_bed_colnames = ["CHROM", "START", "END", "GENE_NAME"]
SNPS_results_colnames = ["CHROM", "POS", "REF", "ALT", "TGT", "GENE_NAME"]


def derive_snps_position(snps_df_with_ref_info: pd.DataFrame) -> pd.DataFrame:
    snps_df_with_ref_info["GENE_POS"] = (
        snps_df_with_ref_info["POS"] - snps_df_with_ref_info["START"]
    )
    return snps_df_with_ref_info


def main(
    ref: str = typer.Option(..., help="Reference file in BED format"),
    snps: str = typer.Option(..., help="File with SNPs information in tab format"),
    output: str = typer.Option(..., help="Filename of the final TSV file"),
):
    ref_df = pd.read_csv(ref, sep="\t", names=REF_bed_colnames)
    snps_df = pd.read_csv(snps, sep="\t", names=SNPS_results_colnames)
    ref_df.drop("CHROM", axis=1, inplace=True)
    snps_df_with_ref_info = snps_df.merge(ref_df, on="GENE_NAME", how="left")
    snps_df_with_snps_position = derive_snps_position(snps_df_with_ref_info)
    snps_df_with_snps_position[[*SNPS_results_colnames, "GENE_POS"]].to_csv(
        output, sep="\t", index=False
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s][%(name)s] -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    typer.run(main)
