Data exploration
================
Geovanny Risco
May 27, 2023

- <a href="#1-import-libraries" id="toc-1-import-libraries">1 Import
  libraries</a>
- <a href="#2-constantconfig-variables"
  id="toc-2-constantconfig-variables">2 Constant/Config variables</a>
- <a href="#3-import-data" id="toc-3-import-data">3 Import data</a>
- <a href="#4-clean-and-prepare-data" id="toc-4-clean-and-prepare-data">4
  Clean and prepare data</a>
  - <a href="#41-args" id="toc-41-args">4.1 ARGs</a>
    - <a href="#411-null-values-detection"
      id="toc-411-null-values-detection">4.1.1 Null values detection</a>
    - <a href="#412-outliers-detection" id="toc-412-outliers-detection">4.1.2
      Outliers detection</a>
  - <a href="#42-feature-engineering" id="toc-42-feature-engineering">4.2
    Feature engineering</a>
  - <a href="#43-snps" id="toc-43-snps">4.3 SNPs</a>
    - <a href="#431-preparation" id="toc-431-preparation">4.3.1
      Preparation</a>
    - <a href="#432-null-values-detection"
      id="toc-432-null-values-detection">4.3.2 Null values detection</a>
    - <a href="#433-outliers-analysis" id="toc-433-outliers-analysis">4.3.3
      Outliers analysis</a>
  - <a href="#44-snps-from-card" id="toc-44-snps-from-card">4.4 SNPs from
    CARD</a>
  - <a href="#45-filtering" id="toc-45-filtering">4.5 Filtering</a>
  - <a href="#46-explorationvisualization"
    id="toc-46-explorationvisualization">4.6 Exploration/Visualization</a>
  - <a href="#47-preparation" id="toc-47-preparation">4.7 Preparation</a>
  - <a href="#48-amr-labels" id="toc-48-amr-labels">4.8 AMR labels</a>
    - <a href="#481-preparation" id="toc-481-preparation">4.8.1
      Preparation</a>
    - <a href="#482-cleaning" id="toc-482-cleaning">4.8.2 Cleaning</a>
    - <a href="#483-exploration" id="toc-483-exploration">4.8.3
      Exploration</a>
- <a href="#5-explore-data" id="toc-5-explore-data">5 Explore data</a>
- <a href="#6-correlation-analysis" id="toc-6-correlation-analysis">6
  Correlation analysis</a>
- <a href="#7-save-data" id="toc-7-save-data">7 Save data</a>

# 1 Import libraries

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 4.1.3

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.2 --
    ## v ggplot2 3.4.1     v purrr   1.0.1
    ## v tibble  3.1.8     v dplyr   1.1.0
    ## v tidyr   1.3.0     v stringr 1.5.0
    ## v readr   2.1.4     v forcats 1.0.0

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tibble' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'readr' was built under R version 4.1.3

    ## Warning: package 'purrr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## Warning: package 'stringr' was built under R version 4.1.3

    ## Warning: package 'forcats' was built under R version 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

# 2 Constant/Config variables

``` r
batch_number <- "_batch3"
MAX_NUMBER_OF_SNPS <- 10
MAX_NULLS_PER_ANTIBIOTIC <- 30 # In percentage
```

# 3 Import data

``` r
# ARGs (Antibiotic Resistance Genes)
args_data_filepath <- paste0("data/results/args_calling/args_table", batch_number, ".csv")
args_data <- read_csv(args_data_filepath, col_types = cols("sample_name" = col_character()))
## Fix/extract name of samples
args_data <- args_data %>%
  mutate(sample_name = str_extract(sample_name, "^\\w+.\\d+"))
args_data
```

    ## # A tibble: 6,306 x 162
    ##    sample_name    NarA  NarB aac(3~1 aac(3~2 aac(6~3 aac(6~4 aadA1 aadA2 aadE-~5
    ##    <chr>         <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl> <dbl> <dbl>   <dbl>
    ##  1 GCA_01263718~     0     0       0       0       0       0     0     0       0
    ##  2 GCA_01263728~     0     0       0       0       0       0     0     0       0
    ##  3 GCA_01263731~     0     0       0       0       0       0     0     0       0
    ##  4 GCA_01263738~     0     0       0       0       0       0     0     0       0
    ##  5 GCA_01263742~     0     0       0       0       0       0     0     0       0
    ##  6 GCA_01263744~     0     0       0       0       0       0     0     0       0
    ##  7 GCA_01263748~     0     0       0       0       0       0     0     0       0
    ##  8 GCA_01263786~     0     0       0       0       0       0     0     0       0
    ##  9 GCA_01264252~     0     0       0       0       0       0     0     0       0
    ## 10 GCA_01264264~     0     0       0       0       0       0     0     0       0
    ## # ... with 6,296 more rows, 152 more variables: `ant(2'')-Ia` <dbl>,
    ## #   `ant(6)-Ia` <dbl>, `aph(2'')-Ia` <dbl>, `aph(2'')-Ic` <dbl>,
    ## #   `aph(2'')-If` <dbl>, `aph(2'')-Ig` <dbl>, `aph(3'')-Ib` <dbl>,
    ## #   `aph(3')-III` <dbl>, `aph(3')-IIa` <dbl>, `aph(3')-Ia` <dbl>,
    ## #   `aph(6)-Ic` <dbl>, `aph(6)-Id` <dbl>, `blaCMY-2` <dbl>, `blaHERA-3` <dbl>,
    ## #   `blaOXA-184` <dbl>, `blaOXA-193` <dbl>, `blaOXA-448` <dbl>,
    ## #   `blaOXA-449` <dbl>, `blaOXA-450` <dbl>, `blaOXA-451` <dbl>, ...

``` r
# SNPs (Single Nucleotide Polymorphisms)
snps_data_filepath <- paste0("data/results/variant_calling/snps_data", batch_number, ".tsv")
snps_data <- read_tsv(snps_data_filepath, col_names = c("chrom", "pos", "ref", "alt", "tgt", "gene_name", "gene_pos", "tax_id", "sample_name"), na = c("."))
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 22136 Columns: 9
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (6): chrom, ref, alt, tgt, gene_name, sample_name
    ## dbl (3): pos, gene_pos, tax_id
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Fix/extract name of samples
snps_data <- snps_data %>%
  mutate(sample_name = str_extract(sample_name, "^\\w+.\\d+"))
snps_data
```

    ## # A tibble: 22,136 x 9
    ##    chrom             pos ref   alt   tgt   gene_name   gene_pos tax_id sample_~1
    ##    <chr>           <dbl> <chr> <chr> <chr> <chr>          <dbl>  <dbl> <chr>    
    ##  1 NZ_KB944666.1 2170478 A     G     A/G   WMS_RS13655      234   1351 1351.853 
    ##  2 NZ_KB944666.1 2170481 C     T     C/T   WMS_RS13655      237   1351 1351.853 
    ##  3 NZ_KB944666.1 2170493 T     C     T/C   WMS_RS13655      249   1351 1351.853 
    ##  4 NZ_KB944666.1 2170505 T     G     T/G   WMS_RS13655      261   1351 1351.853 
    ##  5 NZ_KB944666.1 2170514 G     A     G/A   WMS_RS13655      270   1351 1351.853 
    ##  6 NZ_KB944666.1 2170548 A     G     A/G   WMS_RS13655      304   1351 1351.853 
    ##  7 NZ_KB944666.1 2170556 C     T     C/T   WMS_RS13655      312   1351 1351.853 
    ##  8 NZ_KB944666.1 2170559 T     C     T/C   WMS_RS13655      315   1351 1351.853 
    ##  9 NZ_KB944666.1 2170562 C     G     C/G   WMS_RS13655      318   1351 1351.853 
    ## 10 NZ_KB944666.1 2170571 A     G     A/G   WMS_RS13655      327   1351 1351.853 
    ## # ... with 22,126 more rows, and abbreviated variable name 1: sample_name

``` r
# Results from CARD database
card_data_filepath <- paste0("data/results/card/card_results", "_batch1", ".tsv")
card_data <- read_tsv(card_data_filepath, na = c("n/a"))
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 8120 Columns: 27
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (16): SAMPLE_ID, ORF_ID, Contig, Orientation, Cut_Off, Best_Hit_ARO, Mod...
    ## dbl  (9): TAX_ID, Start, Stop, Pass_Bitscore, Best_Hit_Bitscore, Best_Identi...
    ## lgl  (2): Other_SNPs, Nudged
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Fix/extract name of samples
card_data <- card_data %>%
  mutate(SAMPLE_ID = str_extract(SAMPLE_ID, "^\\w+.\\d+"))
card_data
```

    ## # A tibble: 8,120 x 27
    ##    TAX_ID SAMPLE_ID  ORF_ID Contig  Start   Stop Orien~1 Cut_Off Pass_~2 Best_~3
    ##     <dbl> <chr>      <chr>  <chr>   <dbl>  <dbl> <chr>   <chr>     <dbl>   <dbl>
    ##  1    195 GCA_00528~ AACMV~ AACMV~  40782  41555 +       Perfect     500    516.
    ##  2    195 GCA_00528~ AACMV~ AACMV~   1086   3005 +       Strict      300   1322.
    ##  3    195 GCA_00528~ AACMV~ AACMV~  22802  23575 +       Perfect     500    516.
    ##  4    195 GCA_00528~ AACMV~ AACMV~   1107   3026 +       Strict      300   1322.
    ##  5    195 GCA_00528~ AACMV~ AACMV~ 223686 226253 +       Strict     1200   1496.
    ##  6    195 GCA_00528~ AACMR~ AACMR~  39847  40620 +       Perfect     500    516.
    ##  7    195 GCA_00528~ AACMR~ AACMR~  18355  20274 +       Strict      300   1304.
    ##  8    195 GCA_00528~ AACMR~ AACMR~   3297   5216 +       Strict      300   1304.
    ##  9    195 GCA_00528~ AACMR~ AACMR~   2577   3350 +       Perfect     500    516.
    ## 10    197 GCA_00529~ AACNR~ AACNR~  20452  21225 -       Perfect     500    516.
    ## # ... with 8,110 more rows, 17 more variables: Best_Hit_ARO <chr>,
    ## #   Best_Identities <dbl>, ARO <dbl>, Model_type <chr>,
    ## #   SNPs_in_Best_Hit_ARO <chr>, Other_SNPs <lgl>, `Drug Class` <chr>,
    ## #   `Resistance Mechanism` <chr>, `AMR Gene Family` <chr>, Predicted_DNA <chr>,
    ## #   Predicted_Protein <chr>, CARD_Protein_Sequence <chr>,
    ## #   `Percentage Length of Reference Sequence` <dbl>, ID <chr>, Model_ID <dbl>,
    ## #   Nudged <lgl>, Note <chr>, and abbreviated variable names ...

``` r
# Reference genomes (in BED format)

## Loop through reference genomes folder and read BED files
refs_genomes_dir <- "data/reference_genomes"
refs_genomes_files <- list.files(refs_genomes_dir, pattern = ".bed", full.names = TRUE, recursive = TRUE)
## Read each of the BED files and store them ina dictionary which key is the name of the file and the values a list of genes
refs_genomes <- list()
for (file in refs_genomes_files) {
  tax_id <- str_remove(basename(file), ".bed")
  refs_genomes[tax_id] <- read_tsv(file, col_select = 4)
}
```

    ## Rows: 2727 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): WMS_RS03090
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 1749 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): dcuC
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 1627 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): dnaA
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 4672 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): thrL
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 4493 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): thrL
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# AMR (Antimicrobial Resistance) labels
amr_labels <- read_csv("data/results/data_collection_ncbi/amr_labels.csv", col_types = cols("SampleID" = col_character()))
amr_labels
```

    ## # A tibble: 6,242 x 23
    ##    SampleID     amikacin amoxi~1 ampic~2 cefox~3 cefti~4 ceftr~5 chlor~6 cipro~7
    ##    <chr>           <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 SAMN04256112        0       1       1       1       1       1       0       0
    ##  2 SAMN04256111        0       0       0       0       0       0       0       0
    ##  3 SAMN04256110        0       0       0       0       0       0       0       0
    ##  4 SAMN04256109        0       0       0       0       0       0       0       0
    ##  5 SAMN04256108        0       0       0       0       0       0       0       0
    ##  6 SAMN04256107        0       0       1       0       0       0       0       0
    ##  7 SAMN04256106        0       0       0       0       0       0       0       0
    ##  8 SAMN04256105        0       1       1       1       1       1       0       0
    ##  9 SAMN04256104        0       1       1       1       1       1       0       0
    ## 10 SAMN04256103        0       1       1       1       1       1       0       0
    ## # ... with 6,232 more rows, 14 more variables: gentamicin <dbl>,
    ## #   kanamycin <dbl>, `nalidixic acid` <dbl>, streptomycin <dbl>,
    ## #   sulfisoxazole <dbl>, tetracycline <dbl>,
    ## #   `trimethoprim-sulfamethoxazole` <dbl>, sulfamethoxazole <dbl>,
    ## #   azithromycin <dbl>, meropenem <dbl>, clindamycin <dbl>, erythromycin <dbl>,
    ## #   florfenicol <dbl>, telithromycin <dbl>, and abbreviated variable names
    ## #   1: `amoxicillin-clavulanic acid`, 2: ampicillin, 3: cefoxitin, ...

``` r
# Load NCBI samples metadata
samples_metadata <- read_tsv("data/results/data_collection_ncbi/assembly_accession_ids+tax_ids.txt", col_names = c("biosample_accession", "assembly_accession", "tax_id"))
```

    ## Rows: 6208 Columns: 3
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (2): biosample_accession, assembly_accession
    ## dbl (1): tax_id
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
samples_metadata
```

    ## # A tibble: 6,208 x 3
    ##    biosample_accession assembly_accession tax_id
    ##    <chr>               <chr>               <dbl>
    ##  1 SAMN03842331        GCA_006813885.2     28901
    ##  2 SAMN03988294        GCA_008404065.2     28901
    ##  3 SAMN05596322        GCA_022569775.1     28901
    ##  4 SAMN04563576        GCA_020381005.1       562
    ##  5 SAMN04588431        GCA_008474925.2     28901
    ##  6 SAMN04530411        GCA_008471245.2     28901
    ##  7 SAMN04536994        GCA_008471305.2     28901
    ##  8 SAMN05440588        GCA_008524845.2     28901
    ##  9 SAMN04605174        GCA_008412045.2     28901
    ## 10 SAMN04964189        GCA_008474065.2     28901
    ## # ... with 6,198 more rows

For now and until BV-BRC is active again, we will filter out 1351
samples (which are not in NCBI). In order to add BV-BRC samples, we need
to parse AMR labels information in a different way:

``` r
amr_labels <- amr_labels %>%
  filter(`SampleID` %in% samples_metadata$biosample_accession)
snps_data <- snps_data %>%
  filter(sample_name %in% samples_metadata$assembly_accession)
card_data <- card_data %>%
  filter(SAMPLE_ID %in% samples_metadata$assembly_accession)
args_data <- args_data %>%
  filter(sample_name %in% samples_metadata$assembly_accession)
```

# 4 Clean and prepare data

First of all, we will need to clean and prepare the data in order to
perform the analysis.

## 4.1 ARGs

This table has the following structure:

| sample_name     | GeneA | GeneB | GeneC |   … |
|:----------------|:------|:------|:------|----:|
| GCA_012637185.1 | 0     | 1     | 0     |   … |
| …               | …     | …     | …     |   … |

For each gene, a boolean value is given dependending on whether the gene
is resistance or not.

- 1: resistance gene
- 0: non-resistance gene

### 4.1.1 Null values detection

``` r
# Count number of nulls per sample
args_data %>%
  mutate(nulls = rowSums(is.na(select(., -sample_name)))) %>%
  select(sample_name, nulls) %>%
  arrange(desc(nulls))
```

    ## # A tibble: 6,208 x 2
    ##    sample_name     nulls
    ##    <chr>           <dbl>
    ##  1 GCA_012637185.1     0
    ##  2 GCA_012637285.1     0
    ##  3 GCA_012637315.1     0
    ##  4 GCA_012637385.1     0
    ##  5 GCA_012637425.1     0
    ##  6 GCA_012637445.1     0
    ##  7 GCA_012637485.1     0
    ##  8 GCA_012637865.1     0
    ##  9 GCA_012642525.1     0
    ## 10 GCA_012642645.1     0
    ## # ... with 6,198 more rows

We have no nulls values for ARGs.

### 4.1.2 Outliers detection

``` r
# For each sample, how many resistance genes are present?
args_data %>%
  mutate(count = rowSums(select(., -sample_name))) %>%
  select(sample_name, count) %>%
  arrange(desc(count))
```

    ## # A tibble: 6,208 x 2
    ##    sample_name     count
    ##    <chr>           <dbl>
    ##  1 GCA_007758465.1    23
    ##  2 GCA_008551155.1    21
    ##  3 GCA_008477615.1    20
    ##  4 GCA_008469665.1    19
    ##  5 GCA_007194375.1    19
    ##  6 GCA_007742235.1    18
    ##  7 GCA_008478965.1    18
    ##  8 GCA_008552695.1    18
    ##  9 GCA_008474545.1    18
    ## 10 GCA_007763135.1    18
    ## # ... with 6,198 more rows

``` r
# Mean number of resistance genes per sample
args_data %>%
  mutate(count = rowSums(select(., -sample_name))) %>%
  summarise(mean(count)) %>%
  pull()
```

    ## [1] 3.627416

``` r
# Boxplot summarizing above information
args_data %>%
  mutate(count = rowSums(select(., -sample_name))) %>%
  ggplot(aes(x = "", y = count)) +
  geom_boxplot() +
  labs(x = "", y = "Number of resistance genes", title = "Distribution of resistance genes per sample") +
  theme_bw()
```

![](figures/unnamed-chunk-11-1.png)<!-- -->

## 4.2 Feature engineering

As all of the columns we have in this table is boolean data, there is
not much to do in terms of feature engineering. Hoewever, we have
noticed that there are some ARGs that does not contribute to resistance
in any sample, that is, they are always 0. Since this does not provide
us any information, we will remove them from the table.

``` r
# Check length of args_data columns
original_ncols <- length(colnames(args_data))

# Filter to only columns that has any other value different than 0
args_data <- args_data %>%
  select_if(function(x) any(x != 0))

removed_ncols <- original_ncols - length(colnames(args_data))
```

After the filtering, we have removed 17 columns.

## 4.3 SNPs

### 4.3.1 Preparation

In this case, we will need to perform preparation steps for each sample,
since the table has a different structure.

| chrom         | pos     | ref | alt | tgt | gene_name     | gene_pos | tax_id |     sample_name |
|:--------------|:--------|:----|:----|:----|:--------------|:---------|:-------|----------------:|
| NC_002695.2   | 1250767 | T   | C   | T/C | ECs_1169      | 511.0    | 562    | GCA_012688215.1 |
| NC_002695.2   | 1250712 | C   | A   | C/A | ECs_1169      | 456.0    | 562    | GCA_012688215.1 |
| NZ_CP046317.1 | 218754  | T   | G   | T/G | FOC43_RS01045 | 603      | 195    | GCA_005283725.1 |
| …             | …       | …   | …   | …   | …             | …        | …      |               … |

Here is a brief description of each column:

- chrom: chromosome in which the SNP is located
- pos: position of the SNP in the chromosome
- ref: reference nucleotide (the one that is expected to be found in the
  sample)
- alt: alternative nucleotide (the one that is actually found in the
  sample, that is, the mutation/SNP)
- tgt: translated genotype, or in other words, the ref and alt
  nucleotides in a single string (ref/alt). \#TODO: double check this
- gene_name: name of the gene in which the SNP is located
- gene_pos: position of the SNP in the gene
- tax_id: taxonomic id of the sample
- sample_name: name of the sample

Most of this information has been extracted from the VCF file generated
by the variant calling pipeline. More information about the VCF format
can be found [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

In order to use this data for ML purposes, we need to transform it into
a table with the following structure:

sample_name \| gene_name1/gene1_pos1 \| gene_name1/gene1_pos2 \|
gene_name2/gene2_pos1 \| …

For the previous example, the table would look like this:

| sample_name     | ECs_1169/456 | ECs_1169/1250712 | FOC43_RS01045/603 |   … |
|:----------------|:-------------|:-----------------|:------------------|----:|
| GCA_012688215.1 | C            | NO-SNP           | NULL              |   … |
| GCA_012637185.1 | NO-SNP       | A                | NULL              |   … |
| GCA_005283725.1 | NULL         | NULL             | G                 |   … |

As we can see, each SNP is represented as a column, and each row
represents a sample. For each SNP, we can have the following values:

- NO-SNP: no SNP found in this position of the gene, but the gene does
  belong to the sample
- NULL: the gene does not belong to the sample
- A, C, G, T: the SNP found in this position of the gene

In addition, we will encode these values into numerical values so that
we do limit the number of ML algorithms we can use. But before that,
let’s do some sanity checks.

``` r
# Filter out SNPs with no gene name (which are the same as sample with non gene_pos)
snps_data <- snps_data %>%
  filter(!is.na(gene_name))
# Filter out SNPS which count for >5 in the same gene (this is considered an anomaly, due to reference genome, etc.) #TODO: double check this
snps_data <- snps_data %>%
  group_by(sample_name, gene_name) %>%
  mutate(count = n()) %>%
  filter(count <= MAX_NUMBER_OF_SNPS) %>%
  select(-count)
```

Pivot table with above specifications:

``` r
snps_data_wide <- snps_data %>%
  select(sample_name, gene_name, gene_pos, alt) %>%
  pivot_wider(names_from = c("gene_pos"), values_from = alt, values_fill = "NO-SNP") %>%
  pivot_longer(cols = -c(sample_name, gene_name), names_to = "gene_pos", values_to = "alt") %>%
  pivot_wider(names_from = c("gene_name", "gene_pos"), values_from = alt, values_fill = "NULL", names_sep = "/") %>%
  ungroup()
```

Now we will encode the values into numerical values. We will use the
following encoding:

``` r
snp_2_num <- c(
  "NULL" = -1,
  "NO-SNP" = 0,
  "A" = 1,
  "C" = 2,
  "G" = 3,
  "T" = 4
)
```

``` r
# Mutate all except sample_name
snps_data_wide <- snps_data_wide %>%
  mutate(across(-c(sample_name), ~ as.numeric(snp_2_num[.x])))
```

### 4.3.2 Null values detection

``` r
# Count number of nulls per sample in percentage
snps_data_wide %>%
  mutate(nulls = rowSums(select(., -sample_name) == -1) / ncol(select(., -sample_name)) * 100) %>%
  select(sample_name, nulls) %>%
  arrange(desc(nulls))
```

    ## # A tibble: 114 x 2
    ##    sample_name     nulls
    ##    <chr>           <dbl>
    ##  1 GCA_005287105.1  98.0
    ##  2 GCA_012687565.1  98.0
    ##  3 GCA_005284005.1  98.0
    ##  4 GCA_012708885.1  98.0
    ##  5 GCA_012714385.1  98.0
    ##  6 GCA_012714465.1  98.0
    ##  7 GCA_012708785.1  98.0
    ##  8 GCA_005285885.1  98.0
    ##  9 GCA_005282525.1  98.0
    ## 10 GCA_005285105.1  98.0
    ## # ... with 104 more rows

For most of the samples, there is very little coocurrences in terms of
SNPs.

### 4.3.3 Outliers analysis

``` r
# For each sample, how many SNPs are present?
snps_data_wide %>%
  mutate(count = rowSums(select(., -sample_name) > 0)) %>%
  select(sample_name, count) %>%
  arrange(desc(count))
```

    ## # A tibble: 114 x 2
    ##    sample_name     count
    ##    <chr>           <dbl>
    ##  1 GCA_012686285.1    37
    ##  2 GCA_012688045.1    33
    ##  3 GCA_012688215.1    25
    ##  4 GCA_012717195.1    24
    ##  5 GCA_012686445.1    21
    ##  6 GCA_012735855.1    20
    ##  7 GCA_012714445.1    20
    ##  8 GCA_008524145.1    20
    ##  9 GCA_012642525.1    18
    ## 10 GCA_005289765.1    17
    ## # ... with 104 more rows

``` r
# Mean number of SNPs per sample
snps_data_wide %>%
  mutate(count = rowSums(select(., -sample_name) > 0)) %>%
  summarise(mean(count)) %>%
  pull()
```

    ## [1] 5.850877

``` r
# Boxplot summarizing above information
snps_data_wide %>%
  mutate(count = rowSums(select(., -sample_name) > 0)) %>%
  ggplot(aes(x = "", y = count)) +
  geom_boxplot() +
  labs(x = "", y = "Number of SNPs", title = "Distribution of SNPs per sample") +
  theme_bw()
```

![](figures/unnamed-chunk-20-1.png)<!-- -->

## 4.4 SNPs from CARD

Although [CARD database](https://card.mcmaster.ca/) offers us a large
variety of information about AMR vectors, we will only use the SNPs
information. For more information about the output format, please refer
to the official [documentation](https://github.com/arpcard/rgi#id72).

## 4.5 Filtering

We will be filtering by the following criteria: \* Column `Model_type`
must be either `protein variant model` or `protein overexpression model`
\* They must have a value within the column `SNPs_in_Best_Hit_ARO`.
NOTE: this column can have multiple values separated by commas.

``` r
# Filter by Model_type
card_snps_data <- card_data %>%
  filter(Model_type %in% c("protein variant model", "protein overexpression model"))

# Filter by SNPs_in_Best_Hit_ARO
card_snps_data <- card_snps_data %>%
  filter(!is.na(SNPs_in_Best_Hit_ARO))

# Explode SNPs_in_Best_Hit_ARO
card_snps_data <- card_snps_data %>%
  mutate(SNPs_in_Best_Hit_ARO = strsplit(SNPs_in_Best_Hit_ARO, ",")) %>%
  unnest(SNPs_in_Best_Hit_ARO)
```

## 4.6 Exploration/Visualization

``` r
# Boxplot showing how many SNPs are present in each sample
card_snps_data %>%
  group_by(SAMPLE_ID) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = "", y = count)) +
  geom_boxplot() +
  labs(x = "", y = "Number of SNPs", title = "Distribution of SNPs per sample") +
  theme_bw()
```

![](figures/unnamed-chunk-22-1.png)<!-- -->

## 4.7 Preparation

Now that we have filtered the data, we will need to transform it into a
format compatible for ML algorithms, that is, a table with the features
of interest as columns and the samples as rows. In this case, the
features we are interested in are the SNPs, so we will need to pivot the
table so that each row represents a sample and each column represents a
SNP ID (column `SNPs_in_Best_Hit_ARO`). The value of each cell will be
the number of times that the SNP appears in the sample.

``` r
# Pivot table
card_snps_data_wide <- card_snps_data %>%
  select(SAMPLE_ID, SNPs_in_Best_Hit_ARO) %>%
  group_by(SAMPLE_ID, SNPs_in_Best_Hit_ARO) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = SNPs_in_Best_Hit_ARO, values_from = count, values_fill = 0) %>%
  ungroup()
```

    ## `summarise()` has grouped output by 'SAMPLE_ID'. You can override using the
    ## `.groups` argument.

## 4.8 AMR labels

The structure of this table is as follows:

| SampleID     | Antibiotic1 | Antibiotic2 | Antibiotic3 |   … |
|:-------------|:------------|:------------|:------------|----:|
| SAMN04256112 | 0           | 1           | 0           |   … |
| …            | …           | …           | …           |   … |

Where each column represents an antibiotic and each row represents a
sample. The values of each cell can be:

- 0: the sample is not resistant to the antibiotic
- 1: the sample is resistant to the antibiotic

One sample can be resistant to multiple antibiotics, so we can have
multiple 1s in the same row.

### 4.8.1 Preparation

Adapt data so it has the same sampleIds as ARGS and variant calling
data. AMR labes happens to have the biosamples accession numbers as
sampleIds, so we will need to map them to their corresponding assembly
accession ids.

``` r
# Given information in samples_metadata, replace biosample_accession with assembly_accession
amr_labels <- amr_labels %>%
  left_join(distinct(samples_metadata, biosample_accession, .keep_all = TRUE), by = c("SampleID" = "biosample_accession")) %>%
  select(-c(`SampleID`, tax_id)) %>%
  rename(`SampleID` = assembly_accession) %>%
  select(SampleID, everything())
```

### 4.8.2 Cleaning

We will remove those antibiotics with more than 30% of null values.

``` r
# Count null values per antibiotic (each column) in percentage
nulls_per_antibiotic <- amr_labels %>%
  select(-`SampleID`) %>%
  summarise_all(~ sum(is.na(.x)) / nrow(amr_labels) * 100) %>%
  gather(key = "antibiotic", value = "% of null values") %>%
  arrange(desc(`% of null values`))
knitr::kable(nulls_per_antibiotic)
```

| antibiotic                    | % of null values |
|:------------------------------|-----------------:|
| sulfamethoxazole              |       94.4417593 |
| meropenem                     |       89.5762848 |
| clindamycin                   |       87.6752054 |
| erythromycin                  |       87.6752054 |
| florfenicol                   |       87.6752054 |
| telithromycin                 |       87.6752054 |
| amikacin                      |       54.8251974 |
| azithromycin                  |       45.1909135 |
| kanamycin                     |       37.7315934 |
| ceftiofur                     |       22.7646206 |
| sulfisoxazole                 |       17.8991461 |
| streptomycin                  |       12.3409054 |
| trimethoprim-sulfamethoxazole |       12.3409054 |
| amoxicillin-clavulanic acid   |       12.3247946 |
| ampicillin                    |       12.3247946 |
| cefoxitin                     |       12.3247946 |
| ceftriaxone                   |       12.3247946 |
| chloramphenicol               |       12.3247946 |
| gentamicin                    |        0.0161108 |
| nalidixic acid                |        0.0161108 |
| tetracycline                  |        0.0161108 |
| ciprofloxacin                 |        0.0000000 |

``` r
# Remove antibiotics with more than 30% of null values
antibiotics_to_remove <- nulls_per_antibiotic %>%
  filter(`% of null values` > 30) %>%
  pull(antibiotic)
amr_labels <- amr_labels %>%
  select(-all_of(antibiotics_to_remove))
```

### 4.8.3 Exploration

Count how many samples are resistance to each antibiotic:

``` r
resistant_samples_per_antibiotic <- amr_labels %>%
  select(-`SampleID`) %>%
  summarise_all(~ sum(.x == 1, na.rm = TRUE)) %>%
  gather(key = "antibiotic", value = "resistant samples") %>%
  arrange(desc(`resistant samples`))
knitr::kable(resistant_samples_per_antibiotic)
```

| antibiotic                    | resistant samples |
|:------------------------------|------------------:|
| tetracycline                  |              3343 |
| streptomycin                  |              1818 |
| ampicillin                    |              1652 |
| sulfisoxazole                 |              1539 |
| ceftriaxone                   |               759 |
| amoxicillin-clavulanic acid   |               753 |
| ceftiofur                     |               707 |
| cefoxitin                     |               669 |
| gentamicin                    |               630 |
| nalidixic acid                |               275 |
| chloramphenicol               |               204 |
| ciprofloxacin                 |               185 |
| trimethoprim-sulfamethoxazole |                76 |

``` r
# Bar plot with number of resistant samples per antibiotic
resistant_samples_per_antibiotic %>%
  ggplot(aes(x = antibiotic, y = `resistant samples`)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Antibiotic", y = "Number of resistant samples")
```

![](figures/unnamed-chunk-28-1.png)<!-- -->

# 5 Explore data

Median number of resistant genes per antibiotic:

``` r
# Boxplot with median number of resistant genes per antibiotic
args_data %>%
  mutate(n_args = rowSums(select(., -sample_name))) %>%
  select(sample_name, n_args) %>%
  left_join(amr_labels, by = c("sample_name" = "SampleID")) %>%
  pivot_longer(cols = -c(sample_name, n_args), names_to = "antibiotic", values_to = "resistant") %>%
  filter(resistant == 1) %>%
  ggplot(aes(x = antibiotic, y = n_args, color = antibiotic)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Antibiotic", y = "Number of resistant genes")
```

![](figures/unnamed-chunk-29-1.png)<!-- -->

Median number of SNPs per antibiotic:

``` r
# Boxplot with median number of SNPs per antibiotic
# In this case, SNPS can have multiple values, not only 0 and 1. We will count those with values > 0 as valid SNPs conferring resistance
snps_data %>%
  ungroup() %>%
  select(sample_name) %>%
  group_by(sample_name) %>%
  summarise(n_snps = n()) %>%
  left_join(amr_labels, by = c("sample_name" = "SampleID")) %>%
  pivot_longer(cols = -c(sample_name, n_snps), names_to = "antibiotic", values_to = "resistant") %>%
  filter(resistant == 1) %>%
  ggplot(aes(x = antibiotic, y = n_snps, color = antibiotic)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Antibiotic", y = "Number of SNPs")
```

![](figures/unnamed-chunk-30-1.png)<!-- -->

Median number of CARD SNPs per antibiotic:

``` r
# Boxplot with median number of CARD SNPs per antibiotic
card_snps_data %>%
  group_by(SAMPLE_ID) %>%
  summarise(n_card_snps = n()) %>%
  left_join(amr_labels, by = c("SAMPLE_ID" = "SampleID")) %>%
  pivot_longer(cols = -c(SAMPLE_ID, n_card_snps), names_to = "antibiotic", values_to = "resistant") %>%
  filter(resistant == 1) %>%
  ggplot(aes(x = antibiotic, y = n_card_snps, color = antibiotic)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Antibiotic", y = "Number of CARD SNPs")
```

![](figures/unnamed-chunk-31-1.png)<!-- -->

# 6 Correlation analysis

In the next section we will perform a correlation analysis between the
different variables in the dataset and the multiple antibiotics.

``` r
# Prepate data for correlation analysis: Remove null values and sort data in same order
arranged_amr_labels <- amr_labels %>%
  drop_na() %>%
  arrange(SampleID)
arranged_args_data <- args_data %>%
  filter(sample_name %in% arranged_amr_labels$`SampleID`) %>%
  arrange(sample_name) %>%
  select(-sample_name) %>%
  select_if(function(x) any(x != 0))
arranged_amr_labels <- arranged_amr_labels %>%
  select(-`SampleID`)

# Calculate correlation matrix and its p-values
args_correlation_matrix_coefficients <- matrix(NA, nrow = ncol(arranged_args_data), ncol = ncol(arranged_amr_labels), dimnames = list(colnames(arranged_args_data), colnames(arranged_amr_labels)))
args_correlation_matrix_pvalues <- matrix(NA, nrow = ncol(arranged_args_data), ncol = ncol(arranged_amr_labels), dimnames = list(colnames(arranged_args_data), colnames(arranged_amr_labels)))
for (i in 1:ncol(arranged_args_data)) {
  for (j in 1:ncol(arranged_amr_labels)) {
    args_correlation_matrix_coefficients[i, j] <- cor.test(arranged_args_data[[i]], arranged_amr_labels[[j]])$estimate
    args_correlation_matrix_pvalues[i, j] <- cor.test(arranged_args_data[[i]], arranged_amr_labels[[j]])$p.value
  }
}
#TODO: analyze p-values to see if they are significant
# args_correlation_matrix <- cor(arranged_args_data, arranged_amr_labels) # Alternative method to calculate correlation matrix with coefficients but not p-values

# heatmap of correlation matrix
heatmap(args_correlation_matrix_coefficients,
  xlab = "Antibiotics", ylab = "Antibiotic Resistance Genes (ARGs)",
  main = "Correlation Heatmap",
  col = colorRampPalette(c("blue", "white", "red"))(100),
  key = TRUE,
  key.title = "Correlation Coefficients",
  # Separate axis a bit more
  margins = c(16, 6)
)
```

    ## Warning in plot.window(...): "key" is not a graphical parameter

    ## Warning in plot.window(...): "key.title" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "key" is not a graphical parameter

    ## Warning in plot.xy(xy, type, ...): "key.title" is not a graphical parameter

    ## Warning in title(...): "key" is not a graphical parameter

    ## Warning in title(...): "key.title" is not a graphical parameter

``` r
# Add gradient color legend
legend("topleft",
  legend = c("-1", "0", "1"),
  fill = colorRampPalette(c("blue", "white", "red"))(3),
  title = "Correlation Coefficient",
  cex = 0.8
)
```

![](figures/unnamed-chunk-32-1.png)<!-- -->

# 7 Save data

``` r
snps_data_output_path <- paste0("data/results/variant_calling/snps_data", batch_number, "_cleaned.tsv")
snps_data_wide %>%
  write_tsv(snps_data_output_path)

card_snps_data_output_path <- paste0("data/results/card/card_snps_data", "_batch1", "_cleaned.tsv")
card_snps_data_wide %>%
  write_tsv(card_snps_data_output_path)

args_data_output_path <- paste0("data/results/args_calling/args_data", batch_number, "_cleaned.tsv")
args_data %>%
  write_tsv(args_data_output_path)

amr_labels_output_path <- paste0("data/results/data_collection_ncbi/amr_labels", batch_number, "_cleaned.tsv")
amr_labels %>%
  write_tsv(amr_labels_output_path)
```
