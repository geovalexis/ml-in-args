# Early prediction of Antimicrobial Resistance (AMR) using Machine Learning

## Summary

Using genetic information from resistant bacteria, it is possible to model the characteristics that confer resistance to certain antibiotics. This information could be used to predict new potential resistant genes or sets of genes and take appropriate action. In this project, real data from a series of resistant bacteria identified in the laboratory will be used to build a machine learning model capable of predicting potential resistance genes (ARGS).

## Requirements

* Install [`snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Repository structure

* `data/`: contains input data for the project (not all data is included in the repository due to size limitations)
* `eda/`: exploratory data analysis of the input data
* `src/`: contains the source code of the project
    * `args_calling/`: identifies the ARGs present in the samples by Resfinder v4.3.1
    * `data_retrieval/`: contains various scripts to retrieve the input data
    * `variant_calling/`: variant calling pipeline to get SNPs from the samples
    * `config.yaml`: configuration file for the project (**You might need to modify this file in order to run the project**)

## Future improvements

* Variant calling:
    * Take secondary mutations into account
    * Encode TGT instead of only alternative allele as it provides more information (from what nucleotides the mutation comes from)
    * Use one-hot-encoded instead of label encoding
* General workflow:
    * Use [data-dependant conditional execution](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) to run all workflows instead of custom `run_workflow.py` script
