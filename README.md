
<!-- README.md is generated from README.Rmd. Please edit that file -->

# invivoSyn

<!-- badges: start -->
<!-- badges: end -->

The goal of invivoSyn is to evaluate synergy for in vivo tumor growth
data. Synergy can be calculated based on AUC-based drug effect or linear
mixed model. For effect based efficacy, two reference models can be
selected: HSA (Highest Single Agent) or Bliss (Bliss Independence).
P-values are calculated for both synergistic effect and antagonism
effect.

The package supports **N-drug combinations** (vehicle + N single-agent
arms + 1 combo arm). `read_tv()` accepts explicit role arguments
(`vehicle`, `singles`, `combo`) so the raw CSV can be in any row order,
can contain extra non-participating treatments, and is no longer
restricted to the canonical 4-arm layout.

## Live App

> **<https://hupharm.shinyapps.io/invivoSyn/>** — try the hosted app
> directly in your browser, no local installation required.

## Shiny application

The bundled Shiny application accepts wide or long tumor-volume files,
supports multiple Combination arms with exact comparator mapping, and
runs independent Bliss or HSA TGI/AUC analyses for each Combination.

``` r
invivoSyn::run_invivoSyn_app()
```

For local development, load the package checkout first with
`devtools::load_all()`. For shinyapps.io deployment, install `invivoSyn`
from a reproducible remote source and verify dependencies with
`rsconnect::appDependencies("inst/shiny/invivoSyn")`.

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

# Declare roles explicitly — row order in the CSV no longer matters, extra
# non-participating treatments are filtered out, and the resolved roles travel
# with the tibble for every downstream call.
tv <- read_tv(
  system.file("extdata", "test.csv", package = "invivoSyn"),
  vehicle = "Vehicle",
  singles = c("Rabusertib", "Irinotecan"),
  combo   = "Irinotecan+Rabusertib"
)

# AUC-based efficacy (median AUC ratio, stratified bootstrap)
AUC_lst <- get_mAUCr(tv, ci = 0.9, ci_type = "bca")
AUC_lst$bsAUC_df
#>     Group             Treatment       mAUCr    std.err          lb          ub
#> 1 Group 2            Rabusertib  0.94088183 0.18143161  0.72331064  1.29410019
#> 2 Group 3            Irinotecan  0.01874466 0.04828123 -0.07392685  0.08908393
#> 3 Group 4 Irinotecan+Rabusertib -0.25852501 0.08069597 -0.44118538 -0.15798998

# AUC-based synergy — Bliss and HSA reference models
bliss_synergy_AUC <- AUC_synergy(AUC_lst)
bliss_synergy_AUC
#>          Metric     Value  std.err         lb        ub p.val.Synergy
#> 1            CI 0.5110688 0.272841  0.2148070  1.291781         0.077
#> 2 Synergy_score 2.8224391 5.607127 -0.3361242 19.879188         0.077
#>   p.val.Antagonism
#> 1            0.923
#> 2            0.923
hsa_synergy_AUC   <- AUC_synergy(AUC_lst, method = "HSA", t = 21, ci = 0.9, ci_type = "bca")
hsa_synergy_AUC
#>          Metric     Value    std.err        lb        ub p.val.Synergy
#> 1            CI 0.4302568 0.07467981 0.3201738 0.5647867             0
#> 2 Synergy_score 3.9066779 1.86962557 1.8905558 7.7472475             0
#>   p.val.Antagonism
#> 1                1
#> 2                1
```

### N-drug combinations

`singles` can be any character vector of length ≥ 2 — for a three-drug
combination supply three single-agent names and the combo arm:

``` r
tv <- read_tv("triple_study.csv",
              vehicle = "PBS",
              singles = c("A", "B", "C"),
              combo   = "A+B+C")
AUC_synergy(get_mAUCr(tv))   # Bliss expected survival = prod(S_i)
```

All downstream functions (`AUC_synergy`, `lmm_synergy`, `CombPDX_CI`,
`global_CI_synergy`) iterate over the N singles automatically.
`lmm_synergy` warns when N \> 3, since the full N-way interaction grows
as 2^N − 1.
