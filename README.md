# Early prediction of Antimicrobial Resistance (AMR) using Machine Learning

## Summary

Using genetic information from resistant bacteria, it is possible to model the characteristics that confer resistance to certain antibiotics. This information could be used to predict new potential resistant genes or sets of genes and take appropriate action. In this project, real data from a series of resistant bacteria identified in the laboratory will be used to build a machine learning model capable of predicting potential resistance genes (ARGS).

## Requirements

* Install [`snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Repository structure

* `data/`: contains input data for the project (not all data is included in the repository due to size limitations)
* `eda/`: exploratory data analysis of the input data
* `ml/`: contains all the machine learning related work
* `src/`: contains the source code of the project
    * `card`: identifies ARGs and SNPs present in the samples by [CARD](https://github.com/arpcard/rgi).
    * `data_collection_bvbrc/`: contains various scripts to retrieve the input data from [BV-BRC](https://www.bv-brc.org/) source
    * `data_collection_ncbi/`: contains various scripts to retrieve the input data from [NCBI](https://www.ncbi.nlm.nih.gov/) source
    * `resfinder/`: identifies the ARGs present in the samples by [Resfinder v4.3.1](https://bitbucket.org/genomicepidemiology/resfinder/src/master/src/resfinder/)
    * `config.yaml`: configuration file for the project (**You might need to modify this file in order to run the project**)

## Future improvements

* Obtain more data from [BV-BRC](https://www.bv-brc.org/)
* ARGs calling:
    * Use other ARGs identification tools such as [CARD](https://card.mcmaster.ca/) or [ARG-ANNOT](https://omictools.com/arg-annot-tool)