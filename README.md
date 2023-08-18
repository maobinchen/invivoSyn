
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynergyInVivo

<!-- badges: start -->
<!-- badges: end -->

The goal of SynergyInVivo is to evaluate synergy for in vivo tumor
growth data. Synergy can be calculated based on TGI/AUC based drug
effect or linear mixed model. For effect based efficacy, three reference
models can be selected, which are HSA (Highest Single Agent),
Bliss(Bliss Independence) or RA(Response Addivity).

## Installation

You can install the development version of SynergyInVivo from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("maobinchen/SynergyInVivo")
```

## Example

This is a basic example which shows you how to do synergy calculation

``` r
library(SynergyInVivo)
tv <- read_tv(system.file("extdata", "test.csv", package = "SynergyInVivo"))
TGI_lst <- getTGI(tv,17)
bliss_synergy_TGI <- TGI_synergy(TGI_lst)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
AUC_lst <- get_mAUCr(SNU_81, ci = 0.9, ci_type = "bca")
bliss_synergy_AUC <- AUC_synergy(AUC_lst)
```

<img src="man/figures/README-example-2.png" width="100%" />
