
<!-- README.md is generated from README.Rmd. Please edit that file -->

# invivoSyn

<!-- badges: start -->
<!-- badges: end -->

The goal of invivoSyn is to evaluate synergy for in vivo tumor growth
data. Synergy can be calculated based on TGI/AUC based drug effect or
linear mixed model. For effect based efficacy, three reference models
can be selected, which are HSA (Highest Single Agent), Bliss(Bliss
Independence) or RA(Response Addivity). P-values are calculated for both
synergistic effect and antagonism effect.

The package supports **N-drug combinations** (vehicle + N single-agent
arms + 1 combo arm). `read_tv()` accepts explicit role arguments
(`vehicle`, `singles`, `combo`) so the raw CSV can be in any row order,
can contain extra non-participating treatments, and is no longer
restricted to the canonical 4-arm layout.

## Installation

You can install the development version of invivoSyn from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("maobinchen/invivoSyn")
```

## Example

This is a basic example which shows you how to do synergy calculation

``` r
library(invivoSyn)

# Legacy positional layout (vehicle first, combo last) still works:
tv <- read_tv(system.file("extdata", "test.csv", package = "invivoSyn"))

# Recommended: declare roles explicitly. Row order in the CSV no longer matters,
# extra non-participating treatments are filtered out automatically, and the
# resolved roles travel with the tibble for every downstream call.
tv <- read_tv(
  system.file("extdata", "test.csv", package = "invivoSyn"),
  vehicle = "Vehicle",
  singles = c("Rabusertib", "Irinotecan"),
  combo   = "Irinotecan+Rabusertib"
)

TGI_lst <- getTGI(tv, 17)
bliss_synergy_TGI <- TGI_synergy(TGI_lst)
```

<img src="man/figures/README-example-1.png" alt="" width="100%" />

``` r
TGI_lst_RTV <- getTGI(tv, 17, tv_var = "RTV") # TGI definition from CombPDX paper
bliss_synergy_TGI_RTV <- TGI_synergy(TGI_lst_RTV)
```

<img src="man/figures/README-example-2.png" alt="" width="100%" />

``` r
#global_CI <- global_CI_synergy(tv)
AUC_lst <- get_mAUCr(SNU_81, ci = 0.9, ci_type = "bca")
bliss_synergy_AUC <- AUC_synergy(AUC_lst)
```

### N-drug combinations

`singles` can be any character vector of length ≥ 2 — for a three-drug
combination supply three single-agent names and the combo arm:

``` r
tv <- read_tv("triple_study.csv",
              vehicle = "PBS",
              singles = c("A", "B", "C"),
              combo   = "A+B+C")
TGI_synergy(getTGI(tv, sel_day = 21))   # Bliss expected TGI = 1 - prod(1 - TGI_i)
```

All downstream functions (`TGI_synergy`, `AUC_synergy`, `lmm_synergy`,
`CombPDX_CI`, `global_CI_synergy`) iterate over the N singles
automatically. `lmm_synergy` warns when N \> 3, since the full N-way
interaction grows as 2^N − 1.
