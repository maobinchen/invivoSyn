groups_ui <- function(id) {
  ns <- shiny::NS(id)
  return(bslib::layout_columns(
    col_widths = c(5, 7),
    bslib::card(bslib::card_header("Confirm arm roles"), shiny::uiOutput(ns("roles_ui"))),
    bslib::card(bslib::card_header("Map exact comparators per Combination"), shiny::uiOutput(ns("comparators_ui")))
  ))
}

groups_server <- function(id, tv) {
  shiny::moduleServer(id, function(input, output, session) {
    suggested <- shiny::reactive(invivoSyn:::suggest_arm_roles(unique(tv()$Treatment)))
    output$roles_ui <- shiny::renderUI({
      purrr::map2(suggested()$arm_id, suggested()$Treatment, function(arm, label) {
        selected <- suggested()$role[suggested()$arm_id == arm]
        shiny::selectInput(session$ns(paste0("role_", arm)), label,
          c("Vehicle", "Single agent", "Combination"), selected
        )
      })
    })
    role_map <- shiny::reactive({
      data <- suggested()
      data$role <- purrr::map_chr(data$arm_id, ~ input[[paste0("role_", .x)]] %||% "Single agent")
      return(data)
    })
    output$comparators_ui <- shiny::renderUI({
      roles <- role_map()
      suggested_map <- invivoSyn:::suggest_comparator_map(roles)
      singles <- roles |>
        dplyr::filter(.data$role == "Single agent")
      combos <- roles |>
        dplyr::filter(.data$role == "Combination")
      purrr::map2(combos$arm_id, combos$Treatment, function(arm, label) {
        shiny::selectizeInput(
          session$ns(paste0("comparators_", arm)), label,
          choices = stats::setNames(singles$arm_id, singles$Treatment),
          selected = suggested_map$comparator_arm_id[
            suggested_map$combination_arm_id == arm
          ],
          multiple = TRUE
        )
      })
    })
    comparator_map <- shiny::reactive({
      roles <- role_map()
      combos <- roles$arm_id[roles$role == "Combination"]
      dplyr::bind_rows(purrr::map(combos, function(combo) {
        values <- input[[paste0("comparators_", combo)]] %||% character()
        return(tibble::tibble(combination_arm_id = combo, comparator_arm_id = values))
      }))
    })
    return(list(role_map = role_map, comparator_map = comparator_map))
  })
}
