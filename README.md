
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Temperature, species identity and morphological traits predict carbonate excretion and mineralogy in tropical reef fishes

[![License: MIT + file
LICENSE](https://img.shields.io/badge/License-MIT%20+%20file%20LICENSE-blue.svg)](https://choosealicense.com/licenses/mit/)

The goal of this project is to reproduce all analyses, figures, and
tables of the paper

> Title: Temperature, species identity and morphological traits predict
> carbonate excretion and mineralogy in tropical reef fishes.
>
> Authors: Ghilardi M, Salter MA, Parravicini V, Ferse SCA, Rixen T,
> Wild C, Birkicht M, Perry CT, Berry A, Wilson RW, Mouillot D, Bejarano
> S

## Content

This repository is structured as follow:

-   [:file_folder:
    data/](https://github.com/mattiaghilardi/FishCaCO3Model/tree/master/data):
    contains all data used in the analyses.

-   [:file_folder:
    analysis/](https://github.com/mattiaghilardi/FishCaCO3Model/tree/master/analysis):
    contains R scripts to reproduce the analyses.

-   [:file_folder:
    outputs/](https://github.com/mattiaghilardi/FishCaCO3Model/tree/master/outputs):
    contains all the outputs created during the analyses, including
    figures, tables and models.

-   [:file_folder:
    R/](https://github.com/mattiaghilardi/FishCaCO3Model/tree/master/R):
    contains R functions developed mainly for this project.

-   [:file_folder:
    man/](https://github.com/mattiaghilardi/FishCaCO3Model/tree/master/man):
    contains help files of R functions.

-   [:page_facing_up:
    DESCRIPTION](https://github.com/mattiaghilardi/FishCaCO3Model/blob/master/DESCRIPTION):
    contains project metadata.

-   [:page_facing_up:
    make.R](https://github.com/mattiaghilardi/FishCaCO3Model/blob/master/make.R):
    main R script to run the entire project.

## Instructions

-   Clone this repository (for those not familiar with GitHub, click on
    the green button `Code` on the project main page on GitHub and then
    on `Download ZIP` to download the entire repository, thus unzip it).

-   Open the `.Rproj` file in RStudio or open an R session with working
    directory set to the root of the project.

-   Reproduce the entire project by running:

``` r
source("make.R")
```

## Notes

-   When this project is launched, the appropriate version of `renv`
    should be automatically downloaded and installed into the project
    library. After this has completed, running `source("make.R")` will
    first restore the project library locally on the machine.
-   All required packages and R functions will be loaded.
-   Some steps of the workflow might take time.
-   The project depends on parallel computation and uses 4 cores. It is
    possible to change the number of cores in the `make.R` file under
    `Global options`. On Windows it must be set to 1 and some steps
    might take very long.

## Working environment

    #> R version 4.1.3 (2022-03-10)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 20.04.5 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    #> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices datasets  utils     methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] compiler_4.1.3    magrittr_2.0.1    htmltools_0.5.1.1 tools_4.1.3      
    #>  [5] yaml_2.2.1        stringi_1.7.3     rmarkdown_2.10    knitr_1.33       
    #>  [9] stringr_1.4.0     xfun_0.25         digest_0.6.27     rlang_0.4.11     
    #> [13] renv_0.15.4       evaluate_0.14
