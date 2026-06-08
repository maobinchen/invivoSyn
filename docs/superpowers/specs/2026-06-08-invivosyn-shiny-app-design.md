# invivoSyn Shiny App Design

## Purpose

Create a guided Shiny application for one in vivo combination experiment at a
time. The app lets users upload tumor-volume data, confirm experimental-arm
roles, review data quality, calculate TGI- and AUC-based synergy, and download a
complete report.

The first release supports one vehicle arm, at least two single-agent arms, and
one combination arm. It deliberately excludes LMM and CombPDX analyses so the
same analysis contract applies to two-drug and higher-order combinations.

## Location and Launching

The app will be bundled with the R package under:

```text
inst/shiny/invivoSyn/
```

An exported `run_invivoSyn_app()` function will locate the installed app with
`system.file()` and launch it with `shiny::runApp()`.

During local development, developers load the checkout with
`devtools::load_all()` before launching the app. A shinyapps.io deployment must
install `invivoSyn` from a reproducible remote source, such as GitHub or an
internal CRAN-like repository. The deployment must not depend on an
untraceable local-path package installation.

## Application Structure

The app will borrow the Organoid template's modular Bootstrap 5 structure,
cards, sidebars, interactive tables, plotting patterns, and runtime theme
switching. It will not copy project-specific Organoid analysis logic.

```text
inst/shiny/invivoSyn/
  app.R
  R/
    helpers.R
    mod_upload.R
    mod_groups.R
    mod_review.R
    mod_analysis.R
    mod_report.R
  report_template.Rmd
  www/
```

- `app.R`: navbar, theme selection, module wiring, and session-level state.
- `mod_upload.R`: file upload, sheet selection, format detection, column
  mapping, parsing, and raw preview.
- `mod_groups.R`: role auto-detection and required user confirmation.
- `mod_review.R`: validation results, exclusions, summaries, and growth curves.
- `mod_analysis.R`: explicit-run TGI and AUC synergy calculations and results.
- `mod_report.R`: self-contained HTML report and individual artifact downloads.
- `helpers.R`: pure parsing, validation, role-detection, formatting, plotting,
  and result-interpretation helpers.

## User Workflow

The top-level navigation is:

```text
Upload & Map -> Review -> Analyze -> Report
```

### Upload and Map

The user uploads a CSV, XLS, or XLSX file. The app accepts:

- Wide data: treatment and mouse columns followed by one column per study day.
- Long data: treatment, mouse, day, and tumor-volume columns.

The app auto-detects format and likely columns, then lets the user confirm or
correct the mapping. It proposes experimental roles from treatment names:

- Vehicle
- Two or more single agents
- One Combination

Role detection is advisory. The user must confirm the mapping before continuing.

### Review

The app displays:

- Blocking validation errors and nonblocking warnings
- Parsed and normalized data previews
- Counts for groups, mice, days, and excluded observations
- Interactive tumor-growth curves with selectable `TV`, `RTV`, `DeltaTV`, or
  `logTV`
- The confirmed role mapping

### Analyze

The analysis sidebar contains:

- Metric: TGI or AUC
- Reference model: Bliss, HSA, or RA
- Selected analysis day for TGI
- AUC study end day
- Tumor-volume variable for TGI: `DeltaTV` or `RTV`
- Confidence level
- Bootstrap confidence-interval type
- Bootstrap replicate count
- Explicit **Run Analysis** button

Clicking **Run Analysis** creates an immutable snapshot of normalized data, role
mapping, exclusions, and settings. Later input changes do not silently
recalculate results; they mark existing results as stale until the user runs
the analysis again.

### Report

The user can download:

- A self-contained HTML report
- Each result table as CSV
- Each analysis plot as PNG

## Normalized Data Contract

Both input formats will be converted to a long tibble containing:

```text
Group, Treatment, Mouse, Day, TV, TV0, RTV, DeltaTV, logTV
```

Internal `Group` values will be assigned from confirmed roles, independent of
the uploaded treatment order:

- `Group 1`: Vehicle
- `Group 2` through `Group n-1`: Single agents
- `Group n`: Combination

The app should preserve the original treatment labels for display and reports.

## Validation Rules

Analysis is blocked unless all of these conditions hold:

- Exactly one Vehicle arm is assigned.
- Exactly one Combination arm is assigned.
- At least two Single-agent arms are assigned.
- Every uploaded treatment is assigned exactly one role.
- Treatment-Mouse-Day observations are unique.
- Study days and tumor volumes are numeric.
- Tumor volumes are nonnegative.
- Every mouse has a nonmissing baseline observation.
- The selected TGI day contains Vehicle, Combination, and every Single-agent
  arm.
- Every analysis arm contains enough nonmissing observations for the requested
  bootstrap calculation.

Nonblocking warnings include:

- Missing measurements after baseline
- Unequal group sizes
- Mice observed on different final days
- Zero baseline tumor volume, which makes RTV unavailable
- Very small bootstrap replicate counts
- Role assignments that differ from auto-detection

Errors and warnings must identify the affected groups, mice, days, or columns
and suggest a corrective action.

## Multi-Drug Synergy Calculations

The app supports any number of single-agent arms for TGI and AUC calculations.
The Combination arm is always the last internal group.

### TGI

For single-agent TGI values \(TGI_i\), expected combination TGI is:

- Bliss: \(100 \times [1 - \prod_i(1 - TGI_i / 100)]\)
- HSA: \(\max_i(TGI_i)\)
- RA: \(\sum_i TGI_i\)

The synergy score is:

```text
Observed combination TGI - Expected combination TGI
```

Bootstrap resampling is stratified by group. The app reports the point
estimate, confidence interval, and synergy and antagonism p-values.

### AUC

The app calculates per-mouse normalized AUC and group-level survival-like
effects using the existing package approach. For single-agent effects \(s_i\),
the expected combination effect is:

- Bliss: \(\prod_i s_i\)
- HSA: \(\min_i s_i\)
- RA: \(\sum_i(1 - s_i)\), preserving the package's existing
  response-additivity transformation across all single agents

The implementation must isolate the expected-effect calculation in a tested
helper so model behavior is explicit and can be corrected without changing
Shiny modules.

All calculations and result extraction must use the confirmed arm roles. They
must not rely on fixed row numbers, uploaded treatment order, or an assumption
that the combination is the third nonvehicle result.

### Excluded Methods

LMM, day-specific CombPDX CI, and global CombPDX CI will not appear in the first
release. Their current package implementations assume exactly two single-agent
arms and do not satisfy the approved multi-drug contract.

## Results and Interpretation

The analysis page displays:

- Summary cards for observed effect, expected effect, synergy score,
  confidence interval, and synergy and antagonism p-values
- Tumor-growth curves
- TGI or AUC result tables
- Expected-versus-observed plot
- Bootstrap-distribution plots
- A plain-language result label

Interpretation labels use the synergy-score confidence interval:

- **Synergistic**: interval is entirely above the no-synergy threshold.
- **Antagonistic**: interval is entirely below the no-synergy threshold.
- **Inconclusive**: interval crosses the no-synergy threshold.

The report must state the selected reference model and avoid presenting the
label without the corresponding estimate and uncertainty.

## Report Contents

The self-contained HTML report includes:

- Source filename, generation time, and app/package version
- Uploaded format and confirmed column mapping
- Confirmed experimental-arm roles
- Data-quality errors, warnings, and exclusions
- Analysis settings and snapshot metadata
- Tumor-growth plots
- TGI or AUC tables and plots
- Point estimates, confidence intervals, p-values, and interpretation
- A methods note describing the selected reference model

The report must be reproducible from the analysis snapshot and must not use
newer UI state that has not been rerun.

## Error Handling and Performance

- Parsing and analysis errors are caught and presented as actionable messages.
- Slow bootstraps run only after explicit user action.
- The UI shows progress during parsing, analysis, plotting, and report
  rendering.
- Analysis results are cached within the session by a hash of the analysis
  snapshot and settings.
- Uploaded data and cached results use session-temporary storage and are not
  written into the package directory.
- Report generation is disabled when results are stale or absent.

## Dependencies and Deployment

The app will add its required runtime packages to `DESCRIPTION`, expected to
include `shiny`, `bslib`, `DT`, `plotly`, `readxl`, and report-rendering
dependencies. Existing package analysis dependencies remain authoritative.

Before deployment:

- Run package checks with R 4.5.2 using its explicit executable path.
- Confirm the installed app launches through `run_invivoSyn_app()`.
- Run `rsconnect::appDependencies()` and resolve any dependency with
  `Source: NA`.
- Deploy only from a reproducible package source.

## Verification

Add `testthat` coverage for pure helpers and focused Shiny behavior.

Helper and analysis tests cover:

- Wide CSV parsing
- Wide XLSX parsing
- Long-format parsing
- Column auto-detection and mapping overrides
- Role detection and confirmation
- Every blocking validation rule
- Warning generation
- Two-, three-, and higher-order single-agent TGI calculations
- Two-, three-, and higher-order single-agent AUC calculations
- Bliss, HSA, and RA expected-effect helpers
- Stale-result detection
- Report snapshot integrity

Shiny tests use `shiny::testServer()` for module contracts and reactive state.
A manual end-to-end check uses package sample data and at least one synthetic
three-drug combination dataset.

The final verification sequence is:

1. Run targeted tests.
2. Run the full test suite.
3. Regenerate documentation.
4. Run `devtools::check(args = c("--ignore-vignettes"))`.
5. Launch the app locally with R 4.5.2.
6. Generate and inspect HTML, CSV, and PNG outputs.
7. Run the shinyapps.io dependency preflight.

## Out of Scope

- Multiple Combination arms in one uploaded experiment
- LMM and CombPDX analyses
- Simulation and power analysis
- Persistent user accounts or stored projects
- Cross-experiment comparison
- Automatic analysis recalculation on every input change
