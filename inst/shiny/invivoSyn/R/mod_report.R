report_ui <- function(id) {
  ns <- shiny::NS(id)
  return(bslib::layout_column_wrap(
    width = 1 / 3,
    bslib::card(bslib::card_header("Complete report"), shiny::downloadButton(ns("html"), "Download HTML report")),
    bslib::card(
      bslib::card_header("Efficacy and Synergy results"),
      shiny::selectInput(ns("combo"), "Combination", choices = NULL),
      shiny::downloadButton(ns("efficacy_csv"), "Download efficacy CSV"),
      shiny::downloadButton(ns("synergy_csv"), "Download synergy CSV")
    ),
    bslib::card(bslib::card_header("Tumor growth curve"), shiny::downloadButton(ns("growth_png"), "Download growth PNG")),
    bslib::card(
      bslib::card_header("Synergy plot"),
      shiny::downloadButton(ns("bootstrap_png"), "Download package figure PNG")
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
    combo_detail <- function() {
      snap <- require_fresh()
      shiny::req(input$combo)
      return(snap$result$details[[input$combo]])
    }
    shiny::observe({
      snap <- snapshot()
      if (!is.null(snap)) {
        shiny::updateSelectInput(
          session,
          "combo",
          choices = stats::setNames(
            names(snap$result$details),
            names(snap$result$details)
          )
        )
      }
    })
    output$efficacy_csv <- shiny::downloadHandler(
      filename = function() paste0("invivoSyn_", make.names(input$combo), "_", tolower(combo_detail()$metric), "_efficacy.csv"),
      content = function(file) utils::write.csv(combo_detail()$efficacy, file, row.names = FALSE)
    )
    output$synergy_csv <- shiny::downloadHandler(
      filename = function() paste0("invivoSyn_", make.names(input$combo), "_", tolower(combo_detail()$metric), "_synergy.csv"),
      content = function(file) utils::write.csv(as_display_table(combo_detail()$synergy), file, row.names = FALSE)
    )
    output$growth_png <- shiny::downloadHandler(
      filename = function() "invivoSyn_tumor_growth.png",
      content = function(file) {
        snap <- require_fresh()
        ggplot2::ggsave(file, growth_plot(snap$tv, "TV"), width = 11, height = 7, dpi = 300)
      }
    )
    output$bootstrap_png <- shiny::downloadHandler(
      filename = function() paste0("invivoSyn_", make.names(input$combo), "_", tolower(combo_detail()$metric), "_figure.png"),
      content = function(file) file.copy(combo_detail()$figure, file, overwrite = TRUE)
    )
    output$html <- shiny::downloadHandler(
      filename = function() "invivoSyn_report.html",
      content = function(file) {
        snap <- require_fresh()
        env <- new.env(parent = environment())
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
