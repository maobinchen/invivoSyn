#' Run the invivoSyn Shiny application
#'
#' @param ... Arguments passed to [shiny::runApp()].
#'
#' @return The value returned by [shiny::runApp()].
#' @export
run_invivoSyn_app <- function(...) {
  app_dir <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  if (!nzchar(app_dir)) {
    rlang::abort("The bundled invivoSyn Shiny app could not be found.")
  }
  return(shiny::runApp(app_dir, ...))
}
