review_ui <- function(id) {
  ns <- shiny::NS(id)
  return(shiny::tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::card(
        class = "py-1",
        bslib::card_body(shiny::tags$small("Arms"), shiny::tags$div(class = "fs-5 fw-semibold", shiny::textOutput(ns("n_arms"))))
      ),
      bslib::card(
        class = "py-1",
        bslib::card_body(shiny::tags$small("Combinations"), shiny::tags$div(class = "fs-5 fw-semibold", shiny::textOutput(ns("n_combos"))))
      ),
      bslib::card(
        class = "py-1",
        bslib::card_body(shiny::tags$small("Mice"), shiny::tags$div(class = "fs-5 fw-semibold", shiny::textOutput(ns("n_mice"))))
      ),
      bslib::card(
        class = "py-1",
        bslib::card_body(shiny::tags$small("Study days"), shiny::tags$div(class = "fs-5 fw-semibold", shiny::textOutput(ns("n_days"))))
      )
    ),
    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::layout_column_wrap(
        width = 1,
        heights_equal = "row",
        bslib::card(bslib::card_header("Endpoint Summary"), DT::DTOutput(ns("summary"))),
        bslib::card(
          full_screen = TRUE, bslib::card_header("Tumor-growth curves"),
          shiny::selectInput(ns("y"), "Display", c("TV", "RTV", "DeltaTV", "logTV")),
          plotly::plotlyOutput(ns("growth"))
        )
      ),
      bslib::card(
        bslib::card_header("Validation"),
        DT::DTOutput(ns("issues"))
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
