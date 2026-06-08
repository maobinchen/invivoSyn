# invivoSyn Shiny App Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a bundled, deployable Shiny app that ingests wide or long tumor-volume data, maps multiple Combination arms to exact comparators, runs independent Bliss/HSA TGI and AUC analyses, and downloads reports and artifacts.

**Architecture:** Put testable ingestion, validation, mapping, analysis, ranking, and plotting helpers in package `R/` files. Keep the bundled `inst/shiny/invivoSyn/` app modular and have it call those package helpers. Use immutable analysis snapshots so reports and downloads match the displayed results.

**Tech Stack:** R 4.5.2, invivoSyn, shiny, bslib, DT, plotly, readxl, rmarkdown, testthat.

---

### Task 1: Package setup and launcher

**Files:**
- Modify: `DESCRIPTION`
- Create: `R/run_app.R`
- Create: `tests/testthat.R`
- Create: `tests/testthat/test-run-app.R`

- [ ] Add Shiny runtime packages to `Imports`, test packages to `Suggests`, and testthat edition configuration.
- [ ] Write a failing test that `system.file("shiny/invivoSyn/app.R", package = "invivoSyn")` exists after installation.
- [ ] Add exported `run_invivoSyn_app()` that calls `shiny::runApp(system.file("shiny/invivoSyn", package = "invivoSyn"), ...)`.
- [ ] Run `devtools::document()` and the launcher test with explicit R 4.5.2.
- [ ] Commit package setup.

### Task 2: Ingestion, normalization, and role mapping helpers

**Files:**
- Create: `R/shiny_data.R`
- Create: `tests/testthat/test-shiny-data.R`

- [ ] Write failing tests for long normalization, wide normalization, column-name suggestion, role suggestion, stable arm IDs, and exact comparator maps.
- [ ] Implement `normalize_tv_long()`, `normalize_tv_wide()`, `suggest_tv_columns()`, `suggest_arm_roles()`, and `build_comparator_map()`.
- [ ] Ensure normalized data contains `arm_id`, `Group`, `Treatment`, `Mouse`, `Day`, `TV`, `TV0`, `RTV`, `DeltaTV`, and `logTV`.
- [ ] Run the focused tests and commit.

### Task 3: Validation helpers

**Files:**
- Create: `R/shiny_validation.R`
- Create: `tests/testthat/test-shiny-validation.R`

- [ ] Write failing tests for one Vehicle, multiple Combinations, incomplete/duplicate comparator maps, duplicate observations, missing baseline, invalid values, missing selected-day arms, and warnings.
- [ ] Implement `validate_invivosyn_experiment()` returning `errors`, `warnings`, and `valid`.
- [ ] Make messages identify affected arms, mice, days, or columns.
- [ ] Run focused tests and commit.

### Task 4: Multi-Combination TGI/AUC analysis

**Files:**
- Create: `R/shiny_analysis.R`
- Create: `tests/testthat/test-shiny-analysis.R`

- [ ] Write failing tests for Bliss/HSA expected effects, two- and three-drug combinations, independent multiple-Combination subsets, result ranking, and snapshot hashing.
- [ ] Implement `expected_tgi()`, `expected_auc_effect()`, `analyze_tgi_combination()`, `analyze_auc_combination()`, `analyze_combinations()`, `rank_combination_results()`, and `analysis_snapshot_id()`.
- [ ] Bootstrap within each arm, calculate confidence intervals and synergy/antagonism p-values, and return tidy result and bootstrap tables.
- [ ] Run focused tests and commit.

### Task 5: Bundled modular Shiny UI

**Files:**
- Create: `inst/shiny/invivoSyn/app.R`
- Create: `inst/shiny/invivoSyn/R/helpers.R`
- Create: `inst/shiny/invivoSyn/R/mod_upload.R`
- Create: `inst/shiny/invivoSyn/R/mod_groups.R`
- Create: `inst/shiny/invivoSyn/R/mod_review.R`
- Create: `inst/shiny/invivoSyn/R/mod_analysis.R`
- Create: `inst/shiny/invivoSyn/R/mod_report.R`
- Create: `tests/testthat/test-shiny-modules.R`

- [ ] Build Bootstrap 5 navbar tabs `Upload & Map`, `Review`, `Analyze`, and `Report`, plus runtime theme selection.
- [ ] Implement CSV/XLS/XLSX upload, wide/long detection, column mapping, and previews.
- [ ] Implement role confirmation and exact multi-select comparator mapping per Combination.
- [ ] Implement review cards, validation messages, tables, and interactive growth curves.
- [ ] Implement explicit-run analysis, stale-result detection, ranked summaries, Combination detail selection, tables, and plots.
- [ ] Add focused `shiny::testServer()` tests for mapping, snapshots, and stale-state behavior.
- [ ] Run focused tests and commit.

### Task 6: Reports and downloads

**Files:**
- Create: `inst/shiny/invivoSyn/report_template.Rmd`
- Modify: `inst/shiny/invivoSyn/R/mod_report.R`
- Create: `tests/testthat/test-shiny-report.R`

- [ ] Write a failing snapshot-integrity test.
- [ ] Implement self-contained HTML report rendering from the immutable analysis snapshot.
- [ ] Add CSV downloads for summary and per-Combination results and PNG downloads for growth and analysis plots.
- [ ] Disable report generation when results are missing or stale.
- [ ] Run focused tests and commit.

### Task 7: Documentation and end-to-end verification

**Files:**
- Modify: `README.Rmd`
- Modify: `README.md`
- Modify: `NAMESPACE`
- Create/Modify: generated `man/run_invivoSyn_app.Rd`

- [ ] Document launching the app locally and deploying from a reproducible package source.
- [ ] Run all tests with `"C:\Program Files\R\R-4.5.2\bin\Rscript.exe"`.
- [ ] Run `roxygen2::roxygenise()` and `devtools::check(args = c("--ignore-vignettes"))`.
- [ ] Launch the app locally, verify the main workflow, and generate HTML/CSV/PNG artifacts.
- [ ] Run `rsconnect::appDependencies("inst/shiny/invivoSyn")` when rsconnect is installed and record any deployment blocker.
- [ ] Commit final documentation and generated files.
