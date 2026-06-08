analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  sidebar <- bslib::sidebar(
    shiny::radioButtons(ns("metric"), "Metric", c("TGI", "AUC")),
    shiny::radioButtons(ns("method"), "Reference model", c("Bliss", "HSA"), inline = TRUE),
    shiny::numericInput(ns("selected_day"), "TGI analysis day", 21, min = 0),
    shiny::numericInput(ns("end_day"), "AUC end day", 21, min = 1),
    shiny::radioButtons(ns("tv_var"), "TGI variable", c("DeltaTV", "RTV"), inline = TRUE),
    shiny::selectInput(ns("conf"), "Confidence level", c("90%" = 0.9, "95%" = 0.95), 0.95),
    shiny::numericInput(ns("boot_n"), "Bootstrap replicates", 1000, min = 100, step = 100),
    shiny::actionButton(ns("run"), "Run Analysis", class = "btn-primary w-100"),
    shiny::textOutput(ns("state"))
  )
  return(bslib::layout_sidebar(
    sidebar = sidebar,
    bslib::card(bslib::card_header("Ranked Combination summary"), DT::DTOutput(ns("summary"))),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Combination detail"),
        shiny::selectInput(ns("combo"), "Combination", choices = NULL),
        DT::DTOutput(ns("detail"))
      ),
      bslib::card(full_screen = TRUE, bslib::card_header("Bootstrap distribution"), plotly::plotlyOutput(ns("bootstrap")))
    )
  ))
}

analysis_server <- function(id, tv, role_map, comparator_map, validation) {
  shiny::moduleServer(id, function(input, output, session) {
    settings <- shiny::reactive(list(
      metric = input$metric, method = input$method, selected_day = input$selected_day,
      end_day = input$end_day, tv_var = input$tv_var, conf = as.numeric(input$conf),
      boot_n = as.integer(input$boot_n)
    ))
    current_id <- shiny::reactive(analysis_snapshot_id(tv(), role_map(), comparator_map(), settings()))
    analysis_validation <- shiny::reactive({
      day_value <- if (identical(input$metric, "TGI")) input$selected_day else NULL
      validate_invivosyn_experiment(tv(), role_map(), comparator_map(), day_value)
    })
    snapshot <- shiny::eventReactive(input$run, {
      shiny::validate(shiny::need(validation()$valid, "Resolve review validation errors before analysis."))
      shiny::validate(shiny::need(analysis_validation()$valid, "Selected analysis day is not valid for all required arms."))
      shiny::withProgress(message = "Running Combination analyses", value = 0.2, {
        result <- analyze_combinations(tv(), role_map(), comparator_map(), settings())
        shiny::incProgress(0.7)
        return(list(
          id = current_id(), tv = tv(), role_map = role_map(),
          comparator_map = comparator_map(), settings = settings(), result = result
        ))
      })
    })
    stale <- shiny::reactive(!is.null(snapshot()) && !identical(snapshot()$id, current_id()))
    output$state <- shiny::renderText(if (stale()) "Results are stale. Run analysis again." else "Results match current inputs.")
    output$summary <- DT::renderDT({
      shiny::req(snapshot())
      DT::datatable(snapshot()$result$summary, options = list(scrollX = TRUE))
    })
    shiny::observe({
      req(tv(), role_map(), comparator_map())
      if (!identical(input$metric, "TGI")) {
        return()
      }
      combos <- role_map()$arm_id[role_map()$role == "Combination"]
      vehicles <- role_map()$arm_id[role_map()$role == "Vehicle"]
      if (length(vehicles) != 1L || length(combos) == 0L) {
        return()
      }
      vehicle <- vehicles[[1]]
      candidate_days <- purrr::map_dbl(combos, function(combo) {
        comparators <- comparator_map()$comparator_arm_id[
          comparator_map()$combination_arm_id == combo
        ]
        if (length(comparators) == 0L) {
          return(NA_real_)
        }
        latest_common_day(tv(), c(vehicle, comparators, combo))
      })
      candidate_days <- candidate_days[!is.na(candidate_days)]
      if (length(candidate_days) > 0) {
        shiny::updateNumericInput(session, "selected_day", value = min(candidate_days))
      }
    })
    shiny::observe({
      shiny::req(snapshot())
      result <- snapshot()$result$summary
      shiny::updateSelectInput(
        session,
        "combo",
        choices = stats::setNames(result$combination_treatment, result$combination_treatment)
      )
    })
    detail <- shiny::reactive({
      shiny::req(snapshot(), input$combo)
      dplyr::filter(snapshot()$result$summary, .data$combination_treatment == input$combo)
    })
    output$detail <- DT::renderDT(DT::datatable(detail()))
    output$bootstrap <- plotly::renderPlotly({
      shiny::req(snapshot(), input$combo)
      plotly::ggplotly(bootstrap_plot(snapshot()$result$bootstrap, input$combo))
    })
    return(list(snapshot = snapshot, stale = stale))
  })
}
