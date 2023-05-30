# Early prediction of Antimicrobial Resistance (AMR) using Machine Learning

## Summary

Using genetic information from resistant bacteria, it is possible to model the characteristics that confer resistance to certain antibiotics. This information could be used to predict new potential resistant genes or sets of genes and take appropriate action. In this project, real data from a series of resistant bacteria identified in the laboratory will be used to build a machine learning model capable of predicting potential resistance genes (ARGS).

## Repository structure

* `data/`: contains input data for the project (not all data is included in the repository due to size limitations)
* `eda/`: Exploratory Data Analysis (EDA) of the input data
* `ml/`: contains all the Machine Learning (ML) related work
* `workflows/`: contains ETL workflows to obtain the input data for the project
    * `card`: identifies ARGs and SNPs present in the samples by [CARD](https://github.com/arpcard/rgi).
    * `data_collection_bvbrc/`: contains various scripts to retrieve the input data from [BV-BRC](https://www.bv-brc.org/) source
    * `data_collection_ncbi/`: contains various scripts to retrieve the input data from [NCBI](https://www.ncbi.nlm.nih.gov/) source
    * `resfinder/`: identifies the ARGs present in the samples by [Resfinder v4.3.1](https://bitbucket.org/genomicepidemiology/resfinder/src/master/src/resfinder/)
    * `config.yaml`: configuration file for the workflows