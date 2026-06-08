read_uploaded_table <- function(file_info, sheet = NULL) {
  ext <- tolower(tools::file_ext(file_info$name))
  if (ext == "csv") {
    return(utils::read.csv(file_info$datapath, check.names = FALSE))
  }
  if (ext %in% c("xls", "xlsx")) {
    return(readxl::read_excel(file_info$datapath, sheet = sheet %||% 1))
  }
  rlang::abort("Upload a CSV, XLS, or XLSX file.")
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

format_issue_table <- function(x) {
  if (nrow(x) == 0) return(tibble::tibble(Status = "No issues"))
  return(dplyr::rename(x, Status = .data$type, Message = .data$message))
}

growth_plot <- function(tv, y = "TV") {
  return(ggplot2::ggplot(
    tv,
    ggplot2::aes(
      x = .data$Day, y = .data[[y]], color = .data$Treatment,
      group = interaction(.data$arm_id, .data$Mouse),
      text = paste0("Treatment: ", .data$Treatment, "<br>Mouse: ", .data$Mouse)
    )
  ) +
    ggplot2::geom_line(alpha = 0.65) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = y, color = "Treatment"))
}

bootstrap_plot <- function(bootstrap, combination_treatment) {
  data <- dplyr::filter(bootstrap, .data$combination_treatment == combination_treatment)
  return(ggplot2::ggplot(data, ggplot2::aes(x = .data$score)) +
    ggplot2::geom_histogram(bins = 35, fill = "#4C78A8", alpha = 0.75) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#B22222") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Synergy score", y = "Bootstrap count"))
}
