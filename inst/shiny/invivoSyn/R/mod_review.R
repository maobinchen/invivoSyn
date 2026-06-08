review_ui <- function(id) {
  ns <- shiny::NS(id)
  return(shiny::tagList(
    bslib::layout_column_wrap(
      width = 1 / 4,
      bslib::value_box("Arms", shiny::textOutput(ns("n_arms")), theme = "primary"),
      bslib::value_box("Combinations", shiny::textOutput(ns("n_combos")), theme = "info"),
      bslib::value_box("Mice", shiny::textOutput(ns("n_mice")), theme = "success"),
      bslib::value_box("Study days", shiny::textOutput(ns("n_days")), theme = "warning")
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Validation"), DT::DTOutput(ns("issues"))),
      bslib::card(bslib::card_header("Endpoint Summary"), DT::DTOutput(ns("summary"))),
      bslib::card(
        full_screen = TRUE, bslib::card_header("Tumor-growth curves"),
        shiny::selectInput(ns("y"), "Display", c("TV", "RTV", "DeltaTV", "logTV")),
        plotly::plotlyOutput(ns("growth"))
      )
    )
  ))
}

review_server <- function(id, tv, role_map, comparator_map) {
  shiny::moduleServer(id, function(input, output, session) {
    validation <- shiny::reactive(validate_invivosyn_experiment(
      tv(), role_map(), comparator_map(), NULL
    ))
    output$n_arms <- shiny::renderText(dplyr::n_distinct(tv()$arm_id))
    output$n_combos <- shiny::renderText(sum(role_map()$role == "Combination"))
    output$n_mice <- shiny::renderText(nrow(dplyr::distinct(tv(), .data$arm_id, .data$Mouse)))
    output$n_days <- shiny::renderText(dplyr::n_distinct(tv()$Day))
    output$issues <- DT::renderDT(format_dt_table(dplyr::bind_rows(validation()$errors, validation()$warnings)))
    output$summary <- DT::renderDT({
      format_dt_table(summarize_endpoint_wide(tv(), input$y), options = list(scrollX = TRUE, pageLength = 10))
    })
    output$growth <- plotly::renderPlotly(plotly::ggplotly(growth_plot(tv(), input$y), tooltip = "text"))
    return(validation)
  })
}
