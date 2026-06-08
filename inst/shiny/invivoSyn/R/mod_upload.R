upload_ui <- function(id) {
  ns <- shiny::NS(id)
  return(bslib::layout_columns(
    col_widths = c(4, 8),
    bslib::card(
      bslib::card_header("Upload and column mapping"),
      shiny::fileInput(ns("file"), "Tumor-volume data", accept = c(".csv", ".xls", ".xlsx")),
      shiny::uiOutput(ns("sheet_ui")),
      shiny::radioButtons(ns("format"), "Input format", c("Auto" = "auto", "Wide" = "wide", "Long" = "long")),
      shiny::uiOutput(ns("mapping_ui")),
      shiny::actionButton(ns("parse"), "Parse data", class = "btn-primary"),
      shiny::verbatimTextOutput(ns("status"))
    ),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header("Normalized preview"),
      DT::DTOutput(ns("preview"))
    )
  ))
}

upload_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    raw <- shiny::reactive({
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
      suggestion <- invivoSyn:::suggest_tv_columns(data)
      choices <- names(data)
      shiny::tagList(
        shiny::selectInput(session$ns("treatment"), "Treatment column", choices, suggestion$treatment),
        shiny::selectInput(session$ns("mouse"), "Mouse column", choices, suggestion$mouse),
        shiny::selectInput(session$ns("day"), "Day column (long format)", c("Not used" = "", choices), suggestion$day),
        shiny::selectInput(session$ns("tv"), "Tumor-volume column (long format)", c("Not used" = "", choices), suggestion$tv)
      )
    })
    parsed <- shiny::eventReactive(input$parse, {
      data <- raw()
      fmt <- input$format
      if (fmt == "auto") fmt <- if (nzchar(input$day) && nzchar(input$tv)) "long" else "wide"
      if (fmt == "long") {
        return(invivoSyn:::normalize_tv_long(data, input$treatment, input$mouse, input$day, input$tv))
      }
      return(invivoSyn:::normalize_tv_wide(data, input$treatment, input$mouse))
    })
    output$preview <- DT::renderDT(DT::datatable(parsed(), options = list(scrollX = TRUE)))
    output$status <- shiny::renderText({
      data <- parsed()
      paste(nrow(data), "normalized observations across", dplyr::n_distinct(data$Treatment), "arms.")
    })
    return(list(tv = parsed, filename = shiny::reactive(input$file$name)))
  })
}
