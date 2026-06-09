upload_ui <- function(id) {
  ns <- shiny::NS(id)
  mapping_card <- bslib::card(
    full_screen = TRUE,
    bslib::card_header("Upload and column mapping"),
    shiny::fileInput(ns("file"), "Tumor-volume data", accept = c(".csv", ".xls", ".xlsx")),
    shiny::uiOutput(ns("sheet_ui")),
    shiny::selectInput(ns("example"), "Or load an example dataset", c(
      "Choose…" = "",
      "Combination demo (test.csv)" = "test.csv",
      "SW837" = "SW837.csv",
      "LS_1034" = "LS_1034.csv"
    )),
    shiny::actionButton(ns("load_example"), "Load example", class = "btn-outline-primary mb-2"),
    shiny::radioButtons(ns("format"), "Input format", c("Auto" = "auto", "Wide" = "wide", "Long" = "long")),
    shiny::uiOutput(ns("mapping_ui")),
    shiny::actionButton(ns("parse"), "Parse data", class = "btn-primary"),
    shiny::verbatimTextOutput(ns("status"))
  )
  preview_card <- bslib::card(
    full_screen = TRUE,
    bslib::card_header("Normalized preview"),
    DT::DTOutput(ns("preview"))
  )
  return(split_container(mapping_card, preview_card, sizes = c(4, 8)))
}

upload_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Data source is either an uploaded file or a bundled example dataset.
    source_rv <- shiny::reactiveVal(NULL)
    parse_trigger <- shiny::reactiveVal(0L)

    shiny::observeEvent(input$file, source_rv(list(type = "upload")))
    shiny::observeEvent(input$load_example, {
      shiny::req(nzchar(input$example %||% ""))
      path <- system.file("extdata", input$example, package = "invivoSyn")
      shiny::validate(shiny::need(nzchar(path), "Example dataset not found."))
      source_rv(list(type = "example", path = path, name = input$example))
      parse_trigger(parse_trigger() + 1L) # examples load + auto-parse
    })
    shiny::observeEvent(input$parse, parse_trigger(parse_trigger() + 1L))

    raw <- shiny::reactive({
      src <- source_rv()
      shiny::req(src)
      if (identical(src$type, "example")) {
        return(utils::read.csv(src$path, check.names = FALSE))
      }
      shiny::req(input$file)
      read_uploaded_table(input$file, input$sheet)
    })
    output$sheet_ui <- shiny::renderUI({
      shiny::req(input$file)
      if (tolower(tools::file_ext(input$file$name)) %in% c("xls", "xlsx")) {
        return(shiny::selectInput(session$ns("sheet"), "Worksheet", readxl::excel_sheets(input$file$datapath)))
      }
      return(NULL)
    })
    output$mapping_ui <- shiny::renderUI({
      data <- raw()
      suggestion <- suggest_tv_columns(data)
      choices <- names(data)
      shiny::tagList(
        shiny::selectInput(session$ns("treatment"), "Treatment column", choices, suggestion$treatment),
        shiny::selectInput(session$ns("mouse"), "Mouse column", choices, suggestion$mouse),
        shiny::selectInput(session$ns("day"), "Day column (long format)", c("Not used" = "", choices), suggestion$day),
        shiny::selectInput(session$ns("tv"), "Tumor-volume column (long format)", c("Not used" = "", choices), suggestion$tv)
      )
    })
    parsed <- shiny::eventReactive(parse_trigger(), {
      data <- raw()
      src <- source_rv()
      # Examples auto-parse before the mapping inputs settle, so derive columns
      # from the data itself; uploads use the user-chosen mapping.
      if (identical(src$type, "example")) {
        sg <- suggest_tv_columns(data)
        trt <- sg$treatment
        mse <- sg$mouse
        day <- sg$day
        tvv <- sg$tv
        fmt <- "auto"
      } else {
        shiny::req(input$treatment, input$mouse)
        trt <- input$treatment
        mse <- input$mouse
        day <- input$day
        tvv <- input$tv
        fmt <- input$format %||% "auto"
      }
      # suggest_tv_columns() can return NA for absent day/tv columns (wide data);
      # treat NA/NULL as "not used" so the format test is never NA.
      has_col <- function(x) !is.null(x) && !is.na(x) && nzchar(x)
      if (fmt == "auto") {
        fmt <- if (has_col(day) && has_col(tvv)) "long" else "wide"
      }
      if (fmt == "long") {
        shiny::req(day, tvv)
        return(normalize_tv_long(data, trt, mse, day, tvv))
      }
      return(normalize_tv_wide(data, trt, mse))
    })
    output$preview <- DT::renderDT(format_dt_table(parsed(), options = list(scrollX = TRUE)))
    output$status <- shiny::renderText({
      data <- parsed()
      paste(nrow(data), "normalized observations across", dplyr::n_distinct(data$Treatment), "arms.")
    })
    filename <- shiny::reactive({
      src <- source_rv()
      if (!is.null(src) && identical(src$type, "example")) src$name else input$file$name
    })
    return(list(tv = parsed, filename = filename))
  })
}
