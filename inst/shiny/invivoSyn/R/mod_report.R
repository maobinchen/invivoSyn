report_ui <- function(id) {
  ns <- shiny::NS(id)
  return(bslib::layout_column_wrap(
    width = 1 / 3,
    bslib::card(bslib::card_header("Complete report"), shiny::downloadButton(ns("html"), "Download HTML report")),
    bslib::card(bslib::card_header("Summary table"), shiny::downloadButton(ns("summary_csv"), "Download summary CSV")),
    bslib::card(bslib::card_header("Bootstrap table"), shiny::downloadButton(ns("bootstrap_csv"), "Download bootstrap CSV")),
    bslib::card(bslib::card_header("Growth plot"), shiny::downloadButton(ns("growth_png"), "Download growth PNG")),
    bslib::card(
      bslib::card_header("Combination bootstrap plot"),
      shiny::selectInput(ns("combo"), "Combination", choices = NULL),
      shiny::downloadButton(ns("bootstrap_png"), "Download bootstrap PNG")
    )
  ))
}

report_server <- function(id, snapshot, stale, filename) {
  shiny::moduleServer(id, function(input, output, session) {
    require_fresh <- function() {
      shiny::req(snapshot())
      shiny::validate(shiny::need(!stale(), "Run analysis again before downloading."))
      return(snapshot())
    }
    output$summary_csv <- shiny::downloadHandler(
      filename = function() "invivoSyn_combination_summary.csv",
      content = function(file) utils::write.csv(require_fresh()$result$summary, file, row.names = FALSE)
    )
    output$bootstrap_csv <- shiny::downloadHandler(
      filename = function() "invivoSyn_bootstrap_results.csv",
      content = function(file) utils::write.csv(require_fresh()$result$bootstrap, file, row.names = FALSE)
    )
    shiny::observe({
      snap <- snapshot()
      if (!is.null(snap)) {
        shiny::updateSelectInput(
          session,
          "combo",
          choices = stats::setNames(
            snap$result$summary$combination_treatment,
            snap$result$summary$combination_treatment
          )
        )
      }
    })
    output$growth_png <- shiny::downloadHandler(
      filename = function() "invivoSyn_tumor_growth.png",
      content = function(file) {
        snap <- require_fresh()
        ggplot2::ggsave(file, growth_plot(snap$tv, "TV"), width = 11, height = 7, dpi = 300)
      }
    )
    output$bootstrap_png <- shiny::downloadHandler(
      filename = function() paste0("invivoSyn_", make.names(input$combo), "_bootstrap.png"),
      content = function(file) {
        snap <- require_fresh()
        ggplot2::ggsave(
          file, bootstrap_plot(snap$result$bootstrap, input$combo),
          width = 8, height = 6, dpi = 300
        )
      }
    )
    output$html <- shiny::downloadHandler(
      filename = function() "invivoSyn_report.html",
      content = function(file) {
        snap <- require_fresh()
        env <- new.env(parent = globalenv())
        env$snapshot <- snap
        env$source_filename <- filename()
        rendered <- rmarkdown::render(
          "report_template.Rmd", output_file = basename(file),
          output_dir = dirname(file), envir = env, quiet = TRUE
        )
        if (!identical(normalizePath(rendered), normalizePath(file))) file.copy(rendered, file, overwrite = TRUE)
      }
    )
  })
}
