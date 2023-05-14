Dataset exploration
================
Geovanny Risco
May 14, 2023

- <a href="#1-import-libraries" id="toc-1-import-libraries">1 Import
  libraries</a>
- <a href="#2-read-raw-data" id="toc-2-read-raw-data">2 Read raw data</a>
- <a href="#3-explore-data" id="toc-3-explore-data">3 Explore data</a>
  - <a href="#31-args" id="toc-31-args">3.1 ARGs</a>
  - <a href="#32-snps" id="toc-32-snps">3.2 SNPs</a>
  - <a href="#33-amr-labels" id="toc-33-amr-labels">3.3 AMR labels</a>

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

# 2 Read raw data

``` r
# ARGs (Antibiotic Resistance Genes)
args_data <- read_csv("data/results/args_calling/args_table.csv", col_types = cols("sample_name" = col_character()))
args_data
```

    ## # A tibble: 444 x 66
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
    ## # ... with 434 more rows, 56 more variables: `ant(2'')-Ia` <dbl>,
    ## #   `ant(6)-Ia` <dbl>, `aph(2'')-Ia` <dbl>, `aph(2'')-Ic` <dbl>,
    ## #   `aph(2'')-If` <dbl>, `aph(2'')-Ig` <dbl>, `aph(3'')-Ib` <dbl>,
    ## #   `aph(3')-III` <dbl>, `aph(3')-IIa` <dbl>, `aph(3')-Ia` <dbl>,
    ## #   `aph(6)-Ic` <dbl>, `aph(6)-Id` <dbl>, `blaCMY-2` <dbl>, `blaHERA-3` <dbl>,
    ## #   `blaOXA-184` <dbl>, `blaOXA-193` <dbl>, `blaOXA-448` <dbl>,
    ## #   `blaOXA-449` <dbl>, `blaOXA-450` <dbl>, `blaOXA-451` <dbl>, ...

``` r
# SNPs (Single Nucleotide Polymorphisms)
snps_data <- read_tsv("data/results/variant_calling/snps_data.tsv", col_names = c("chrom", "pos", "ref", "alt", "tgt", "gene_name", "gene_pos", "tax_id", "sample_name"), na = c("."))
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 13892 Columns: 9
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (6): chrom, ref, alt, tgt, gene_name, sample_name
    ## dbl (3): pos, gene_pos, tax_id
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
snps_data
```

    ## # A tibble: 13,892 x 9
    ##    chrom             pos ref   alt   tgt   gene_name   gene_pos tax_id sample_~1
    ##    <chr>           <dbl> <chr> <chr> <chr> <chr>          <dbl>  <dbl> <chr>    
    ##  1 NZ_KB944666.1  289637 T     C     T/C   <NA>              NA   1351 1351.819 
    ##  2 NZ_KB944666.1  289643 T     C     T/C   <NA>              NA   1351 1351.819 
    ##  3 NZ_KB944666.1  289665 G     C     G/C   <NA>              NA   1351 1351.819 
    ##  4 NZ_KB944666.1  289816 T     A     T/A   <NA>              NA   1351 1351.819 
    ##  5 NZ_KB944666.1  289853 T     A     T/A   <NA>              NA   1351 1351.819 
    ##  6 NZ_KB944666.1  289887 C     T     C/T   <NA>              NA   1351 1351.819 
    ##  7 NZ_KB944666.1  289920 T     C     T/C   <NA>              NA   1351 1351.819 
    ##  8 NZ_KB944666.1 2172187 A     G     A/G   WMS_RS13660       27   1351 1351.819 
    ##  9 NZ_KB944666.1 2172223 T     G     T/G   WMS_RS13660       63   1351 1351.819 
    ## 10 NZ_KB944666.1 2172255 A     G     A/G   WMS_RS13660       95   1351 1351.819 
    ## # ... with 13,882 more rows, and abbreviated variable name 1: sample_name

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

# 3 Explore data

## 3.1 ARGs

``` r
# TODO: plot number of resistance genes per antibiotic
```

## 3.2 SNPs

``` r
# TODO
```

## 3.3 AMR labels

Count how many samples are resistance to each antibiotic:

``` r
amr_labels %>%
  select(-`SampleID`) %>%
  replace(is.na(.), 0) %>% # NOTE: if no result for a certain antibiotic, it is considered as not resistant
  summarise_all(~ sum(.x == 1)) %>%
  gather(key = "antibiotic", value = "resistant samples") %>%
  arrange(desc(`resistant samples`))
```

    ## # A tibble: 22 x 2
    ##    antibiotic                  `resistant samples`
    ##    <chr>                                     <int>
    ##  1 tetracycline                               3364
    ##  2 streptomycin                               1827
    ##  3 ampicillin                                 1666
    ##  4 sulfisoxazole                              1546
    ##  5 ceftriaxone                                 766
    ##  6 amoxicillin-clavulanic acid                 757
    ##  7 ceftiofur                                   711
    ##  8 cefoxitin                                   673
    ##  9 gentamicin                                  632
    ## 10 kanamycin                                   479
    ## # ... with 12 more rows
