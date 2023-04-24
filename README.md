# Machine Learning (ML) in Antibiotic Resistance Genes (ARGs)

## Summary

From the genetic information of resistant bacteria, it is possible to model the characteristics that confer resistance to certain antibiotics. This information could be used to predict new potential resistant genes (or set of genes) and act accordingly. In this project, real data from a series of resistant bacteria identified in the laboratory will be used to build a Machine Learning model capable of predicting potential resistance genes (ARGs).

## Requirements

* Install conda environment from `environment.yml` file:

```
conda env create -f environment.yml
```

* Install other linux dependencies:
    - samtools v1.17 and bcftools v1.17: http://www.htslib.org/download/
    - vcftools v0.1.16: https://github.com/vcftools/vcftools