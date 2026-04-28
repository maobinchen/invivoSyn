# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`invivoSyn` is an R package for calculating and visualizing drug combination synergy from in vivo tumor volume data. It supports three synergy reference models (Bliss Independence, HSA, Response Additivity) and multiple efficacy metrics (TGI, AUC, Linear Mixed Model, CombPDX CI).

## Common Commands

```r
# Regenerate documentation from roxygen2 tags
roxygen2::roxygenise()

# Check package (vignettes excluded per project config)
devtools::check(args = c("--ignore-vignettes"))

# Build package
devtools::build(args = c("--no-build-vignettes"))

# Load package interactively during development
devtools::load_all()

# Run a specific example / manual test
tv <- invivoSyn::read_tv(system.file("extdata", "test.csv", package = "invivoSyn"))
```

There are no automated unit tests (`tests/` directory is absent). Validation is done manually via the vignette (`vignettes/my-vignette.Rmd`) and sample datasets in `data/` and `inst/extdata/`.

## Architecture

### Data flow

```
CSV (wide format)
  └─ read_tv()          → tv tibble (wide: one column per day)
       └─ expand_tv()   → long tibble with RTV, DeltaTV, logTV, TV0 columns
            ├─ getTGI() / get_AUC() / get_mAUCr()   → efficacy scalars
            └─ *_synergy() functions                 → synergy score + CI
```

### Source files

| File | Responsibility |
|------|---------------|
| `R/read_tv.R` | `read_tv()`, `expand_tv()` — data ingestion & reshape |
| `R/efficacy.R` | `getTGI()`, `get_AUC()`, `get_mAUCr()` — per-metric efficacy |
| `R/synergy.R` | `TGI_synergy()`, `AUC_synergy()`, `lmm_synergy()`, `CombPDX_CI()`, `global_CI_synergy()` — core synergy calculations |
| `R/plot.R` | `plot_tumor_growth_curve()`, `plot_group_by_day()`, `plot_density()` |
| `R/power.R` | `power_calc()`, `sim_power()` — Monte Carlo power analysis |
| `R/simulation.R` | `simu_TV()` — synthetic tumor volume data generator |
| `R/utils.R` | `theme_Publication()`, `getCI()`, color scales |

### Synergy calculation pattern

All `*_synergy()` functions follow the same contract:
- Accept a `tv` tibble (output of `read_tv()`) plus group-name arguments (`combo`, `drug1`, `drug2`, `vehicle`)
- Accept a `method` argument: `"Bliss"`, `"HSA"`, or `"RA"`
- Return a named list with the synergy score and bootstrap confidence interval
- Bootstrap CIs use `boot::boot()` internally; parallel bootstraps use `furrr::future_map()`

### Key design conventions

- Column naming in `expand_tv()` output: `Group`, `Subject`, `Day`, `TV`, `TV0`, `RTV`, `DeltaTV`, `logTV`
- Group labels must match exactly across all `*_synergy()` calls (case-sensitive strings)
- `CombPDX_CI()` and `global_CI_synergy()` implement the CombPDX delta-method approach (day-specific and global, respectively) and require the `nlme` mixed-model fit
