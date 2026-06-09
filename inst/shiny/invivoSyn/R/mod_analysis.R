analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  sidebar <- bslib::sidebar(
    shiny::radioButtons(ns("metric"), "Metric", c("TGI", "AUC")),
    shiny::radioButtons(ns("method"), "Reference model", c("Bliss", "HSA"), inline = TRUE),
    shiny::selectInput(ns("selected_day"), "TGI analysis day", choices = NULL),
    shiny::selectInput(ns("end_day"), "AUC end day", choices = NULL),
    shiny::numericInput(ns("auc_t"), "AUC survival-estimation day", 21, min = 1),
    shiny::radioButtons(ns("tv_var"), "TGI variable", c("DeltaTV", "RTV"), inline = TRUE),
    shiny::selectInput(ns("conf"), "Confidence level", c("90%" = 0.9, "95%" = 0.95), 0.95),
    shiny::numericInput(ns("boot_n"), "Bootstrap replicates", 1000, min = 100, step = 100),
    shiny::actionButton(ns("run"), "Run Analysis", class = "btn-primary w-100"),
    shiny::textOutput(ns("state"))
  )
  results_card <- bslib::card(
    full_screen = TRUE,
    bslib::card_header(shiny::textOutput(ns("results_header"))),
    shiny::selectInput(ns("combo"), "Combination", choices = NULL),
    bslib::navset_card_tab(
      bslib::nav_panel("Efficacy", htmltools::div(style = "height:100%; overflow:auto;", DT::DTOutput(ns("efficacy")))),
      bslib::nav_panel("Synergy", htmltools::div(style = "height:100%; overflow:auto;", DT::DTOutput(ns("synergy"))))
    )
  )
  figure_card <- bslib::card(
    full_screen = TRUE,
    bslib::card_header("Synergy plot"),
    shiny::imageOutput(ns("bootstrap"), height = "100%")
  )
  return(bslib::layout_sidebar(
    sidebar = sidebar,
    split_container(
      results_card, figure_card,
      direction = "vertical", sizes = c(1, 1), height = "calc(100vh - 7rem)"
    )
  ))
}

analysis_server <- function(id, tv, role_map, comparator_map, validation) {
  shiny::moduleServer(id, function(input, output, session) {
    valid_days <- shiny::reactive({
      shiny::req(tv())
      sort(unique(stats::na.omit(tv()$Day)))
    })
    settings <- shiny::reactive(list(
      metric = input$metric,
      method = input$method,
      selected_day = as.numeric(input$selected_day),
      end_day = if (identical(input$end_day, "__all__")) NA_real_ else as.numeric(input$end_day),
      auc_t = as.numeric(input$auc_t),
      tv_var = input$tv_var,
      conf = as.numeric(input$conf),
      boot_n = as.integer(input$boot_n)
    ))
    current_id <- shiny::reactive(analysis_snapshot_id(tv(), role_map(), comparator_map(), settings()))
    analysis_validation <- shiny::reactive({
      validate_invivosyn_experiment(tv(), role_map(), comparator_map(), input$selected_day)
    })
    snapshot <- shiny::eventReactive(input$run, {
      shiny::validate(shiny::need(validation()$valid, "Resolve review validation errors before analysis."))
      shiny::validate(shiny::need(analysis_validation()$valid, "Selected analysis day is not valid for all required arms."))
      shiny::withProgress(message = "Running Combination analyses", value = 0.2, {
        result <- analyze_combinations(tv(), role_map(), comparator_map(), settings())
        shiny::incProgress(0.7)
        return(list(
          id = current_id(),
          tv = tv(),
          role_map = role_map(),
          comparator_map = comparator_map(),
          settings = settings(),
          result = result
        ))
      })
    })
    stale <- shiny::reactive(!is.null(snapshot()) && !identical(snapshot()$id, current_id()))
    output$state <- shiny::renderText(if (stale()) "Results are stale. Run analysis again." else "Results match current inputs.")
    output$results_header <- shiny::renderText({
      metric <- if (!is.null(snapshot())) snapshot()$result$metric else input$metric
      paste0("Synergy calculation results (", metric, ")")
    })
    shiny::observe({
      days <- valid_days()
      shiny::updateSelectInput(
        session,
        "selected_day",
        choices = stats::setNames(as.character(days), as.character(days)),
        selected = as.character(max(days))
      )
      shiny::updateSelectInput(
        session,
        "end_day",
        choices = c("All data" = "__all__", stats::setNames(as.character(days), as.character(days))),
        selected = "__all__"
      )
    })
    shiny::observe({
      shiny::req(tv(), role_map())
      # Default the TGI analysis day to the last day on which every arm still has
      # at least 80% of its baseline animals measured.
      coverage_day <- latest_coverage_day(tv(), role_map()$arm_id, min_prop = 0.8)
      selected <- if (is.na(coverage_day)) max(valid_days()) else coverage_day
      shiny::updateSelectInput(session, "selected_day", selected = as.character(selected))
    })
    shiny::observe({
      shiny::req(snapshot())
      shiny::updateSelectInput(
        session,
        "combo",
        choices = stats::setNames(
          names(snapshot()$result$details),
          names(snapshot()$result$details)
        )
      )
    })
    combo_detail <- shiny::reactive({
      shiny::req(snapshot(), input$combo)
      snapshot()$result$details[[input$combo]]
    })
    output$efficacy <- DT::renderDT({
      shiny::req(combo_detail())
      format_dt_table(combo_detail()$efficacy, options = list(scrollX = TRUE))
    })
    output$synergy <- DT::renderDT({
      shiny::req(combo_detail())
      format_dt_table(as_display_table(combo_detail()$synergy), options = list(scrollX = TRUE))
    })
    output$bootstrap <- shiny::renderImage({
      shiny::req(combo_detail())
      list(src = combo_detail()$figure, contentType = "image/png",
           alt = paste(combo_detail()$metric, "package figure"))
    }, deleteFile = FALSE)
    return(list(snapshot = snapshot, stale = stale))
  })
}
