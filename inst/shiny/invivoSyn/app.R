suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(DT)
  library(plotly)
})

options(shiny.maxRequestSize = 500 * 1024^2)

APP_HELPERS <- c(
  "suggest_tv_columns", "normalize_tv_long", "normalize_tv_wide",
  "suggest_arm_roles", "suggest_comparator_map",
  "validate_invivosyn_experiment", "analysis_snapshot_id",
  "analyze_combinations", "build_analysis_tv", "build_combo_role_args",
  "latest_common_day", "latest_coverage_day", "lookup_treatment", "lookup_arm_id"
)

PACKAGE_APIS <- c(
  "set_roles", "get_roles", "getTGI", "get_mAUCr",
  "TGI_synergy", "AUC_synergy"
)

app_env <- environment()
source_files <- c("shiny_data.R", "shiny_validation.R", "shiny_analysis.R")
source_roots <- c(
  file.path("..", "..", "..", "R"),
  file.path("..", "..", "R")
)
has_source_files <- vapply(
  source_roots,
  function(root) all(file.exists(file.path(root, source_files))),
  logical(1)
)
source_root <- source_roots[has_source_files][1]

if (!is.na(source_root)) {
  purrr::walk(
    source_files,
    ~ sys.source(file.path(source_root, .x), envir = app_env)
  )
}

if (!requireNamespace("invivoSyn", quietly = TRUE)) {
  stop("Install invivoSyn or load the invivoSyn package before running the app.")
}

purrr::walk(c(APP_HELPERS, PACKAGE_APIS), function(helper) {
  if (!exists(helper, envir = app_env, inherits = FALSE)) {
    assign(helper, getFromNamespace(helper, "invivoSyn"), envir = app_env)
  }
})

purrr::walk(
  list.files("R", pattern = "\\.R$", full.names = TRUE),
  ~ sys.source(.x, envir = app_env)
)

app_theme <- bslib::bs_theme(
  version = 5,
  bg = "#FAF8F3",
  fg = "#1F2421",
  primary = "#0F6E5C",
  secondary = "#6B7280",
  success = "#2E7D5B",
  info = "#0F6E5C",
  warning = "#B8860B",
  danger = "#9A3B3B",
  base_font = bslib::font_google("Asap", local = FALSE),
  heading_font = bslib::font_google("Fraunces", local = FALSE),
  code_font = bslib::font_google("IBM Plex Mono", local = FALSE),
  "border-radius" = "0.5rem",
  "navbar-bg" = "#FFFFFF"
)

ui <- bslib::page_navbar(
  id = "nav",
  title = shiny::tagList(shiny::icon("chart-line"), shiny::span("invivoSyn", class = "app-wordmark")),
  theme = app_theme,
  fillable = FALSE,
  header = shiny::tags$head(
    shiny::tags$style(shiny::HTML(
      ".app-wordmark { font-family: var(--bs-heading-font-family, 'Fraunces', Georgia, serif); font-weight: 600; letter-spacing: .005em; }
       .navbar { padding-top: 0.75rem; padding-bottom: 0.75rem; }
       .navbar-brand { font-size: 1.5rem; }
       .navbar-nav .nav-link { font-size: 1.05rem; padding-top: 0.6rem; padding-bottom: 0.6rem; }
       .bslib-card > .card-header { font-family: var(--bs-heading-font-family, 'Fraunces', Georgia, serif); font-weight: 600; font-size: 1.1rem; letter-spacing: .01em; padding-top: .55rem; padding-bottom: .55rem; }
       .stat-caption { text-transform: uppercase; letter-spacing: .08em; font-size: .7rem; color: var(--bs-secondary); }
       .stat-value { font-family: var(--bs-code-font-family, 'IBM Plex Mono', monospace); font-size: 1.6rem; line-height: 1; color: var(--bs-primary); }
       .nav-proceed { border-top: 1px solid var(--bs-border-color); margin-top: .75rem; padding-top: .75rem; }
       .nav-proceed .btn-outline-secondary.disabled { border-style: dashed; opacity: .65; }
       .resizable-frame { border: 1px solid var(--bs-border-color); border-radius: .375rem; padding: .25rem; }
       .resizable-frame > .shiny-plot-output,
       .resizable-frame .plotly.html-widget,
       .resizable-frame .shiny-image-output,
       .resizable-frame .dataTables_wrapper { height: 100% !important; width: 100% !important; }
       .split-container { display: flex; width: 100%; }
       .split-horizontal { flex-direction: row; }
       .split-vertical { flex-direction: column; }
       .split-pane { min-width: 80px; min-height: 80px; }
       .split-pane > .bslib-card { height: 100%; margin-bottom: 0; }
       .split-pane > .split-container { height: 100%; }
       .split-pane > .resizable-frame { height: 100%; }
       .split-pane .card-body { overflow: auto; }
       .split-pane .shiny-image-output { height: 100%; }
       .shiny-image-output img { max-width: 100%; max-height: 100%; width: auto; height: auto; object-fit: contain; display: block; margin: 0 auto; }
       .split-gutter { flex: 0 0 10px; background: var(--bs-border-color); border-radius: 4px; }
       .split-horizontal > .split-gutter { cursor: col-resize; margin: 0 .25rem; }
       .split-vertical > .split-gutter { cursor: row-resize; margin: .25rem 0; }
       .split-gutter:hover { background: var(--bs-primary); }"
    )),
    shiny::tags$script(shiny::HTML(
      "(function () {
         var fire = function () { window.dispatchEvent(new Event('resize')); };
         var t; var debounced = function () { clearTimeout(t); t = setTimeout(fire, 150); };
         var ro = new ResizeObserver(debounced);
         var attach = function (el) { if (el.dataset.roAttached) return; el.dataset.roAttached = '1'; ro.observe(el); };
         var scan = function () { document.querySelectorAll('.resizable-frame').forEach(attach); };
         document.addEventListener('DOMContentLoaded', scan);
         new MutationObserver(scan).observe(document.documentElement, { childList: true, subtree: true });
       })();"
    )),
    shiny::tags$script(shiny::HTML(
      "(function () {
         document.addEventListener('mousedown', function (e) {
           var g = e.target;
           if (!g.classList || !g.classList.contains('split-gutter')) return;
           var container = g.parentElement;
           var horizontal = container.classList.contains('split-horizontal');
           var prev = g.previousElementSibling, next = g.nextElementSibling;
           if (!prev || !next) return;
           e.preventDefault();
           var start = horizontal ? e.clientX : e.clientY;
           var pr = prev.getBoundingClientRect(), nr = next.getBoundingClientRect();
           var prevSize = horizontal ? pr.width : pr.height;
           var nextSize = horizontal ? nr.width : nr.height;
           document.body.style.userSelect = 'none';
           var rt;
           function onMove(ev) {
             var delta = (horizontal ? ev.clientX : ev.clientY) - start;
             var np = Math.max(80, prevSize + delta);
             var nn = Math.max(80, nextSize - delta);
             prev.style.flex = '0 0 ' + np + 'px';
             next.style.flex = '0 0 ' + nn + 'px';
             clearTimeout(rt); rt = setTimeout(function () { window.dispatchEvent(new Event('resize')); }, 60);
           }
           function onUp() {
             document.removeEventListener('mousemove', onMove);
             document.removeEventListener('mouseup', onUp);
             document.body.style.userSelect = '';
             window.dispatchEvent(new Event('resize'));
           }
           document.addEventListener('mousemove', onMove);
           document.addEventListener('mouseup', onUp);
         });
       })();"
    ))
  ),
  bslib::nav_panel(
    "Upload & Map", icon = shiny::icon("upload"),
    upload_ui("upload"), groups_ui("groups"),
    shiny::uiOutput("nav_to_review")
  ),
  bslib::nav_panel(
    "Review", icon = shiny::icon("clipboard-check"),
    review_ui("review"),
    shiny::uiOutput("nav_to_analyze")
  ),
  bslib::nav_panel(
    "Analyze", icon = shiny::icon("flask"),
    analysis_ui("analysis"),
    shiny::uiOutput("nav_to_report")
  ),
  bslib::nav_panel("Report", icon = shiny::icon("file-arrow-down"), report_ui("report")),
  bslib::nav_spacer(),
  bslib::nav_panel(
    "Tutorial",
    icon = shiny::icon("circle-question"),
    shiny::tags$iframe(
      src = "tutorial.html",
      style = "width:100%; height:calc(100vh - 56px); border:none; display:block;"
    )
  )
)

server <- function(input, output, session) {
  uploaded <- upload_server("upload")
  mapped <- groups_server("groups", uploaded$tv)
  validation <- review_server("review", uploaded$tv, mapped$role_map, mapped$comparator_map)
  analysis <- analysis_server("analysis", uploaded$tv, mapped$role_map, mapped$comparator_map, validation)
  report_server("report", analysis$snapshot, analysis$stale, uploaded$filename)

  # ---- Guided stage navigation -------------------------------------------
  # Readiness: data parsed + validation clean gates the move past mapping and
  # review; a fresh analysis snapshot gates the move to reporting. Guarded with
  # tryCatch because reading the not-yet-parsed `tv()` eventReactive raises a
  # silent req() error.
  data_validated <- shiny::reactive({
    out <- FALSE
    tryCatch(
      if (shiny::isTruthy(uploaded$tv())) out <- isTRUE(validation()$valid),
      error = function(e) NULL
    )
    return(out)
  })
  ready_review <- shiny::reactive(data_validated())
  ready_analyze <- shiny::reactive(data_validated())
  ready_report <- shiny::reactive(!is.null(analysis$snapshot()) && !isTRUE(analysis$stale()))

  output$nav_to_review <- shiny::renderUI(next_button("to_review", "Proceed to Review", ready_review()))
  output$nav_to_analyze <- shiny::renderUI(next_button("to_analyze", "Proceed to Analyze", ready_analyze()))
  output$nav_to_report <- shiny::renderUI(next_button("to_report", "Proceed to Report", ready_report()))

  shiny::observeEvent(input$to_review, bslib::nav_select("nav", "Review", session = session))
  shiny::observeEvent(input$to_analyze, bslib::nav_select("nav", "Analyze", session = session))
  shiny::observeEvent(input$to_report, bslib::nav_select("nav", "Report", session = session))

  # One-time toast each time a stage transitions from not-ready to ready.
  notify_on_ready <- function(ready_reactive, message) {
    previous <- shiny::reactiveVal(FALSE)
    shiny::observe({
      ready <- isTRUE(ready_reactive())
      if (ready && !isTRUE(previous())) {
        shiny::showNotification(message, type = "message", duration = 5)
      }
      previous(ready)
    })
  }
  notify_on_ready(ready_review, "Data validated — ready for Review.")
  notify_on_ready(ready_report, "Analysis complete — ready for Report.")

  return(invisible(NULL))
}

shiny::shinyApp(ui, server)
