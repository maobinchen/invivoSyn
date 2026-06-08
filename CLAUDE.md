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

# Built-in sample datasets (lazy-loaded from data/)
data(LS_1034); data(SNU_81)
TGI_lst <- getTGI(LS_1034, 17)
TGI_synergy(TGI_lst, method = "Bliss")
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
| `R/zzz.R` | Package `.onLoad` / namespace hooks |

### Synergy calculation pattern

Synergy function input contract:
- `TGI_synergy()` / `AUC_synergy()` accept an **efficacy list** produced by `getTGI()` / `get_AUC()` / `get_mAUCr()` — not the raw `tv` tibble.
- `lmm_synergy()`, `CombPDX_CI()`, `global_CI_synergy()` accept the `tv` tibble directly (they fit their own `nlme` mixed model internally).
- All accept a `method` argument: `"Bliss"`, `"HSA"`, or `"RA"`.
- Return a list containing the synergy score, bootstrap CI, and p-values for both **synergy and antagonism** effects.
- Bootstrap CIs use `boot::boot()`; parallel bootstraps use `furrr::future_map()` (set a plan via `future::plan(multisession)`).

### Key design conventions

- Column naming in `expand_tv()` output: `Group`, `Treatment`, `Mouse`, `Day`, `TV`, `TV0`, `RTV`, `DeltaTV`, `logTV`
- `read_tv()` accepts explicit role arguments (`vehicle`, `singles`, `combo`) and attaches a `roles` attribute consumed by every downstream function via the `get_roles()` helper. Legacy positional 4-group CSVs still work without those arguments.
- The pipeline is N-drug-generic: `singles` may be any character vector of length ≥ 2 (vehicle + N singles + combo = N+2 groups). `CombPDX_CI` / `bs_global_CI` / `bs_AUC_synergy` apply the closed-form Bliss/HSA/RA formulas generalized to N singles; `lmm_synergy` builds the formula `Day:(d1*d2*...*dN)` programmatically (warns when N > 3 due to the 2^N − 1 interaction explosion).
- `CombPDX_CI()` and `global_CI_synergy()` implement the CombPDX delta-method approach (day-specific and global, respectively).
