# invivoSyn Shiny App Design

## Purpose

Create a guided Shiny application for one in vivo experiment at a time. The app
lets users upload tumor-volume data, confirm experimental-arm roles and exact
comparators for each Combination arm, review data quality, calculate TGI- and
AUC-based synergy, and download a complete report.

The first release supports one Vehicle arm, a shared pool of Single-agent arms,
and one or more Combination arms. Combination arms may represent different
component sets, doses, or schedules. It deliberately excludes LMM and CombPDX
analyses so the same analysis contract applies to two-drug and higher-order
combinations.

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
- `mod_groups.R`: role auto-detection, required user confirmation, and exact
  Single-agent comparator mapping for each Combination.
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
- One or more Combinations

Role and component detection is advisory. For every Combination, the app
auto-suggests likely components from treatment names, then requires the user to
map and confirm its exact Single-agent comparator arms. This exact mapping also
handles dose- and schedule-matched comparator selection. Comparator arms may be
reused across Combinations when explicitly mapped.

### Review

The app displays:

- Blocking validation errors and nonblocking warnings
- Parsed and normalized data previews
- Counts for groups, mice, days, and excluded observations
- Interactive tumor-growth curves with selectable `TV`, `RTV`, `DeltaTV`, or
  `logTV`
- The confirmed role and Combination-to-comparator mappings

### Analyze

The analysis sidebar contains:

- Metric: TGI or AUC
- Reference model: Bliss or HSA
- Selected analysis day for TGI
- AUC study end day
- Tumor-volume variable for TGI: `DeltaTV` or `RTV`
- Confidence level
- Bootstrap confidence-interval type
- Bootstrap replicate count
- Explicit **Run Analysis** button

Clicking **Run Analysis** creates an immutable snapshot of normalized data, role
and comparator mappings, exclusions, and settings. Each Combination is analyzed
independently using only the Vehicle, its mapped Single-agent comparators, and
that Combination. Later input changes do not silently recalculate results; they
mark existing results as stale until the user runs the analysis again.

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

The normalized dataset will preserve stable internal arm identifiers and the
original treatment labels. A separate confirmed mapping table will contain:

```text
arm_id, Treatment, role
```

where `role` is `Vehicle`, `Single agent`, or `Combination`. A second mapping
table will contain one row per exact comparator relationship:

```text
combination_arm_id, comparator_arm_id
```

Calculations must derive temporary per-Combination analysis subsets and must not
depend on uploaded treatment order or globally fixed `Group` numbers.

## Validation Rules

Analysis is blocked unless all of these conditions hold:

- Exactly one Vehicle arm is assigned.
- At least one Combination arm is assigned.
- At least two Single-agent arms are assigned.
- Every uploaded treatment is assigned exactly one role.
- Every Combination is mapped to at least two distinct Single-agent comparator
  arms.
- Every mapped comparator exists and has the Single-agent role.
- Every Combination-to-comparator mapping is confirmed by the user.
- Treatment-Mouse-Day observations are unique.
- Study days and tumor volumes are numeric.
- Tumor volumes are nonnegative.
- Every mouse has a nonmissing baseline observation.
- For each Combination, the selected TGI day contains the Vehicle, that
  Combination, and all its mapped Single-agent comparator arms.
- Every arm in each per-Combination analysis subset contains enough nonmissing
  observations for the requested bootstrap calculation.

Nonblocking warnings include:

- Missing measurements after baseline
- Unequal group sizes
- Mice observed on different final days
- Zero baseline tumor volume, which makes RTV unavailable
- Very small bootstrap replicate counts
- Role assignments that differ from auto-detection
- Comparator mappings that differ from auto-detection
- Explicit reuse of a Single-agent comparator across multiple Combinations

Errors and warnings must identify the affected groups, mice, days, or columns
and suggest a corrective action.

## Multi-Drug Synergy Calculations

The app supports one or more Combination arms. Each Combination is analyzed
independently against exactly its confirmed Single-agent comparator arms and
the shared Vehicle. Each per-Combination analysis supports any number of mapped
Single-agent comparators.

### TGI

For single-agent TGI values \(TGI_i\), expected combination TGI is:

- Bliss: \(100 \times [1 - \prod_i(1 - TGI_i / 100)]\)
- HSA: \(\max_i(TGI_i)\)

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

The implementation must isolate the expected-effect calculation in a tested
helper so model behavior is explicit and can be corrected without changing
Shiny modules.

All calculations and result extraction must use the confirmed arm roles and
Combination-to-comparator mappings. They must not rely on fixed row numbers,
uploaded treatment order, or an assumption about a Combination's position.

### Excluded Methods

LMM, day-specific CombPDX CI, and global CombPDX CI will not appear in the first
release. Their current package implementations assume exactly two single-agent
arms and do not satisfy the approved multi-drug contract.

Response Additivity (RA) will also not appear in the first release. The app
supports only Bliss and HSA reference models.

## Results and Interpretation

The analysis page displays:

- A cross-Combination summary table ranked separately within each metric and
  reference model
- Summary cards for observed effect, expected effect, synergy score,
  confidence interval, and synergy and antagonism p-values for the selected
  Combination
- Tumor-growth curves
- Per-Combination TGI or AUC result tables
- Per-Combination expected-versus-observed plots
- Per-Combination bootstrap-distribution plots
- A plain-language result label

The user can select a Combination to inspect its detailed results. The app
ranks Combinations only within the same metric and reference model. It must not
create a single ranking that mixes TGI with AUC or Bliss with HSA. Within each
metric and reference-model pair, Combinations are ranked by synergy score from
highest to lowest. The summary displays confidence intervals and p-values beside
the rank so ranking is not presented as statistical certainty.

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
- Confirmed experimental-arm roles and every Combination-to-comparator mapping
- Data-quality errors, warnings, and exclusions
- Analysis settings and snapshot metadata
- Tumor-growth plots
- Cross-Combination summary tables
- Per-Combination TGI or AUC tables and plots
- Per-Combination point estimates, confidence intervals, p-values, and
  interpretation
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
- Combination component suggestion and exact comparator confirmation
- Every blocking validation rule
- Warning generation
- Two-, three-, and higher-order single-agent TGI calculations
- Two-, three-, and higher-order single-agent AUC calculations
- Multiple Combination arms with different component sets
- Multiple Combination arms with dose- or schedule-specific comparators
- Explicit comparator reuse across Combinations
- Per-Combination analysis isolation
- Cross-Combination ranking within metric and reference model
- Bliss and HSA expected-effect helpers
- Stale-result detection
- Report snapshot integrity

Shiny tests use `shiny::testServer()` for module contracts and reactive state.
A manual end-to-end check uses package sample data and a synthetic experiment
containing multiple Combinations, including different component sets and
dose- or schedule-specific comparator mappings.

The final verification sequence is:

1. Run targeted tests.
2. Run the full test suite.
3. Regenerate documentation.
4. Run `devtools::check(args = c("--ignore-vignettes"))`.
5. Launch the app locally with R 4.5.2.
6. Generate and inspect HTML, CSV, and PNG outputs.
7. Run the shinyapps.io dependency preflight.

## Out of Scope

- LMM and CombPDX analyses
- Simulation and power analysis
- Persistent user accounts or stored projects
- Cross-experiment comparison
- Automatic analysis recalculation on every input change
