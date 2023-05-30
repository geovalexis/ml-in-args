# ETL workflows

## Summary

This folder contains multiple scripts to download and run the necessary tools to obtain the input data for the project.
The scripts are managed and run by [`snakemake`](https://snakemake.readthedocs.io/en/stable/), a workflow manager that 
allows to run reproducible and scalable data analyses. It manages the dependencies between the different scripts and provides
a way to run them in parallell and in the correct order.

## Requirements

* Install [`snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Usage

Before running the workflows, make sure to edit the `config.yaml` file to set the correct paths and config variables. The following variables are required:
* `biosamples_accessions`: CSV file with the biosamples accessions to download from NCBI
* `samples_dir`: directory where the samples will be downloaded

Each of the folders contains a `Snakefile` that can be run with `snakemake` with the following command:

```bash
snakemake --use-conda --cores <number_of_cores>
```

## Future improvements

* Obtain and automate data collection from [BV-BRC](https://www.bv-brc.org/)
* Merge ARGs data from CARD and Resfinder