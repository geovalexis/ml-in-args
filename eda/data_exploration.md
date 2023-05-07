Dataset exploration
================
Geovanny Risco
May 04, 2023

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
args_data <- read_csv("src/args_calling/results/args_table.csv", col_types = cols("sample_name" = col_character()))
head(args_data)
```

    ## # A tibble: 6 x 26
    ##   sample_name  NarA  NarB aac(6'~1 ant(6~2 aph(2~3 aph(2~4 aph(3~5   cat cat(p~6
    ##   <chr>       <dbl> <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl> <dbl>   <dbl>
    ## 1 1351.846        0     0        0       0       0       0       0     0       0
    ## 2 1351.780        0     0        0       0       0       0       0     0       0
    ## 3 1351.812        1     1        1       1       0       0       1     0       0
    ## 4 1351.809        0     0        0       0       0       0       0     0       0
    ## 5 1351.798        0     0        1       1       0       0       1     0       0
    ## 6 1351.791        0     0        0       0       0       0       0     0       0
    ## # ... with 16 more variables: dfrD <dbl>, dfrG <dbl>, `erm(B)` <dbl>,
    ## #   gyrA <dbl>, `lnu(A)` <dbl>, `lnu(B)` <dbl>, `lnu(G)` <dbl>, `lsa(A)` <dbl>,
    ## #   `lsa(E)` <dbl>, parC <dbl>, str <dbl>, `tet(L)` <dbl>, `tet(M)` <dbl>,
    ## #   `tet(O)` <dbl>, `tet(S)` <dbl>, `vat(E)` <dbl>, and abbreviated variable
    ## #   names 1: `aac(6')-aph(2'')`, 2: `ant(6)-Ia`, 3: `aph(2'')-Ia`,
    ## #   4: `aph(2'')-If`, 5: `aph(3')-III`, 6: `cat(pC221)`

``` r
# SNPs (Single Nucleotide Polymorphisms)
snps_data <- read_csv("src/variant_calling/results/snps_table.csv", col_types = cols("sample_name" = col_character()))
head(snps_data)
```

    ## # A tibble: 6 x 315
    ##   sample_name `2170548` 217057~1 21721~2 21721~3 21722~4 21722~5 21722~6 21722~7
    ##   <chr>           <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 1351.846            0        0       0       0       0       0       0       0
    ## 2 1351.780            2        2       4       1       2       1       1       2
    ## 3 1351.812            0        0       0       0       0       0       0       0
    ## 4 1351.809            0        0       0       0       2       1       1       2
    ## 5 1351.798            0        0       0       0       0       0       0       0
    ## 6 1351.791            2        2       0       0       2       1       1       2
    ## # ... with 306 more variables: `2172283` <dbl>, `2172286` <dbl>,
    ## #   `2172295` <dbl>, `2172301` <dbl>, `2172303` <dbl>, `2172307` <dbl>,
    ## #   `2172322` <dbl>, `2172348` <dbl>, `1532357` <dbl>, `1532362` <dbl>,
    ## #   `1532366` <dbl>, `2172176` <dbl>, `2172223` <dbl>, `2172259` <dbl>,
    ## #   `2172312` <dbl>, `2172351` <dbl>, `1532372` <dbl>, `2368677` <dbl>,
    ## #   `2368687` <dbl>, `2368688` <dbl>, `2368692` <dbl>, `2368695` <dbl>,
    ## #   `2368716` <dbl>, `3463` <dbl>, `3481` <dbl>, `3484` <dbl>, ...

``` r
# AMR (Antimicrobial Resistance) labels
amr_labels <- read_csv("src/data_retrieval/results/amr_labels.csv", col_types = cols("Genome ID" = col_character()))
head(amr_labels)
```

    ## # A tibble: 6 x 16
    ##   `Genome ID` chloramp~1 cipro~2 dapto~3 eryth~4 genta~5 kanam~6 linez~7 nitro~8
    ##   <chr>            <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 1351.777             0       0       0       0       0       1       0       0
    ## 2 1351.778             0       0       0       1       1       1       0       0
    ## 3 1351.779             0       0       0       0       0       0       0       0
    ## 4 1351.78              0       0       0       0       0       0       0       0
    ## 5 1351.781             0       0       0       1       0       0       0       0
    ## 6 1351.782             0       0       0       1       0       0       0       0
    ## # ... with 7 more variables: penicillin <dbl>,
    ## #   `quinupristin/dalfopristin` <dbl>, streptomycin <dbl>, tetracycline <dbl>,
    ## #   tigecycline <dbl>, tylosin <dbl>, vancomycin <dbl>, and abbreviated
    ## #   variable names 1: chloramphenicol, 2: ciprofloxacin, 3: daptomycin,
    ## #   4: erythromycin, 5: gentamicin, 6: kanamycin, 7: linezolid,
    ## #   8: nitrofurantoin

# 3 Explore data

## 3.1 ARGs

``` r
#TODO: plot number of resistance genes per antibiotic
```

## 3.2 SNPs

``` r
#TODO
```

## 3.3 AMR labels

Count how many samples are resistance to each antibiotic:

``` r
amr_labels %>%
  select(-`Genome ID`) %>%
  summarise_all(~ sum(.x == 1)) %>%
  gather(key = "antibiotic", value = "resistant samples") %>%
  arrange(desc(`resistant samples`))
```

    ## # A tibble: 15 x 2
    ##    antibiotic                `resistant samples`
    ##    <chr>                                   <int>
    ##  1 tetracycline                               62
    ##  2 tylosin                                    55
    ##  3 erythromycin                               52
    ##  4 kanamycin                                  48
    ##  5 gentamicin                                 30
    ##  6 streptomycin                               28
    ##  7 chloramphenicol                            19
    ##  8 ciprofloxacin                               1
    ##  9 daptomycin                                  1
    ## 10 nitrofurantoin                              1
    ## 11 quinupristin/dalfopristin                   1
    ## 12 linezolid                                   0
    ## 13 penicillin                                  0
    ## 14 tigecycline                                 0
    ## 15 vancomycin                                  0
