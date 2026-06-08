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

format_3 <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.3f", round(x, 3)))
}

format_dt_table <- function(data, ...) {
  out <- DT::datatable(data, ...)
  numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  if (length(numeric_cols) > 0) {
    out <- DT::formatRound(out, columns = numeric_cols, digits = 3)
  }
  return(out)
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
      text = paste0(
        "Treatment: ", .data$Treatment,
        "<br>Mouse: ", .data$Mouse,
        "<br>Day: ", .data$Day,
        "<br>", y, ": ", format_3(.data[[y]])
      )
    )
  ) +
    ggplot2::geom_line(alpha = 0.65) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = y, color = "Treatment"))
}

summarize_endpoint_wide <- function(tv, y = "TV") {
  stats <- tv |>
    dplyr::summarise(
      mean = mean(.data[[y]], na.rm = TRUE),
      sd = stats::sd(.data[[y]], na.rm = TRUE),
      n = sum(!is.na(.data[[y]])),
      .by = c(Group, Treatment, Day)
    ) |>
    dplyr::mutate(
      cell = paste0(
        format_3(.data$mean),
        "+/-",
        ifelse(is.na(.data$sd), "NA", format_3(.data$sd)),
        "(",
        .data$n,
        ")"
      ),
      Day = as.character(.data$Day)
    ) |>
    dplyr::select("Group", "Treatment", "Day", "cell")
  out <- tidyr::pivot_wider(stats, names_from = "Day", values_from = "cell")
  return(as.data.frame(out))
}

bootstrap_plot <- function(bootstrap, combination_treatment) {
  data <- dplyr::filter(bootstrap, .data$combination_treatment == combination_treatment)
  return(ggplot2::ggplot(data, ggplot2::aes(x = .data$score)) +
    ggplot2::geom_histogram(bins = 35, fill = "#4C78A8", alpha = 0.75) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#B22222") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Synergy score", y = "Bootstrap count"))
}
