# ETL workflows

## Summary

This folder contains multiple scripts to download and run the necessary tools to obtain the input data for the project.
The scripts are managed and run by [`snakemake`](https://snakemake.readthedocs.io/en/stable/), a workflow manager that 
allows to run reproducible and scalable data analyses. It manages the dependencies between the different scripts and provides
a way to run them in parallell and in the correct order.

## Requirements

* Install [`snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Usage

Before running the workflows, make sure to edit the `config.yaml` file to set the correct paths and config variables.
The required variables depend on the workflow, but they are all documented in the `config.yaml` file.
Each of the folders contains a `Snakefile` that can be run with `snakemake` with the following command:

```bash
snakemake --use-conda --cores <number_of_cores>
```

Some workflows may support running in the Azure cloud. To do so, you will need to correctly configure the credentials. 
A bash script is provided for that end.

## Future improvements

* Merge ARGs data from CARD and Resfinder