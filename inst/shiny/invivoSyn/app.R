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
  "analyze_combinations"
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
} else {
  if (!requireNamespace("invivoSyn", quietly = TRUE)) {
    stop("Install invivoSyn or run the app from the invivoSyn package checkout.")
  }
  purrr::walk(APP_HELPERS, function(helper) {
    assign(helper, getFromNamespace(helper, "invivoSyn"), envir = app_env)
  })
}

purrr::walk(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

ui <- bslib::page_navbar(
  title = shiny::tagList(shiny::icon("chart-line"), "invivoSyn"),
  theme = bslib::bs_theme(version = 5, bootswatch = "cosmo"),
  fillable = TRUE,
  bslib::nav_panel("Upload & Map", shiny::icon("upload"), upload_ui("upload"), groups_ui("groups")),
  bslib::nav_panel("Review", shiny::icon("clipboard-check"), review_ui("review")),
  bslib::nav_panel("Analyze", shiny::icon("flask"), analysis_ui("analysis")),
  bslib::nav_panel("Report", shiny::icon("file-arrow-down"), report_ui("report")),
  bslib::nav_spacer(),
  bslib::nav_item(shiny::selectInput("theme", NULL, c("Cosmo" = "cosmo", "Flatly" = "flatly", "Darkly" = "darkly"), "cosmo"))
)

server <- function(input, output, session) {
  shiny::observeEvent(input$theme, {
    session$setCurrentTheme(bslib::bs_theme(version = 5, bootswatch = input$theme))
  }, ignoreInit = TRUE)
  uploaded <- upload_server("upload")
  mapped <- groups_server("groups", uploaded$tv)
  selected_day <- shiny::reactive(max(uploaded$tv()$Day, na.rm = TRUE))
  validation <- review_server("review", uploaded$tv, mapped$role_map, mapped$comparator_map, selected_day)
  analysis <- analysis_server("analysis", uploaded$tv, mapped$role_map, mapped$comparator_map, validation)
  report_server("report", analysis$snapshot, analysis$stale, uploaded$filename)
  return(invisible(NULL))
}

shiny::shinyApp(ui, server)
