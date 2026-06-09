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

# Flexbox split with draggable gutters between panes (handled by the splitter JS in
# app.R). `sizes` are flex-grow weights for the initial split; `direction` is
# "horizontal" (drag width) or "vertical" (drag height; supply `height`).
split_container <- function(..., direction = "horizontal", sizes = NULL, height = NULL) {
  panes <- list(...)
  n <- length(panes)
  if (is.null(sizes)) sizes <- rep(1, n)
  kids <- list()
  for (i in seq_len(n)) {
    kids[[length(kids) + 1]] <- htmltools::div(
      class = "split-pane",
      style = sprintf("flex: %s 1 0; overflow: auto;", sizes[[i]]),
      panes[[i]]
    )
    if (i < n) {
      kids[[length(kids) + 1]] <- htmltools::div(class = "split-gutter")
    }
  }
  htmltools::tags$div(
    class = paste("split-container", paste0("split-", direction)),
    style = if (!is.null(height)) sprintf("height: %s;", height) else NULL,
    kids
  )
}

# Wrap content in a user-resizable block. Kept out of bslib's fill flow (place it
# inside card_body(fillable = FALSE)) so its explicit height holds and the native
# CSS resize grip works.
resizable_frame <- function(..., height = "440px") {
  htmltools::div(
    class = "resizable-frame",
    style = sprintf(
      "resize: both; overflow: auto; height: %s; min-height: 180px; width: 100%%;",
      height
    ),
    ...
  )
}

# Stage-navigation button: a full-width footer action. Solid primary when `ready`;
# otherwise a dashed, dimmed outline (reads as a locked next step, not a dead slab)
# with a helper line. Non-clickable via the native `disabled` attribute (no shinyjs).
next_button <- function(id, label, ready) {
  btn <- shiny::actionButton(
    id, label,
    class = if (isTRUE(ready)) "btn-primary btn-lg w-100" else "btn-outline-secondary btn-lg w-100 disabled",
    icon = shiny::icon("arrow-right")
  )
  hint <- NULL
  if (!isTRUE(ready)) {
    btn <- htmltools::tagAppendAttributes(btn, disabled = "disabled", `aria-disabled` = "true")
    hint <- shiny::tags$div(class = "text-muted small text-center mt-1", "Complete this step to continue")
  }
  return(htmltools::div(class = "nav-proceed", btn, hint))
}

# Arm-role colour language for tumour-growth curves: vehicle/control reads as a
# neutral grey; every other arm cycles a curated accent palette led by the app's
# teal. Heuristic name match mirrors suggest_arm_roles()/set_roles().
arm_palette <- function(treatments) {
  trt <- unique(as.character(treatments))
  low <- tolower(trt)
  accent <- c("#0F6E5C", "#9A3B3B", "#B8860B", "#3B5BA5", "#6D4C9F", "#1B7F8C", "#C2603B")
  pal <- character(length(trt))
  names(pal) <- trt
  ai <- 1L
  for (i in seq_along(trt)) {
    if (grepl("vehicle|control|placebo", low[i])) {
      pal[i] <- "#9AA0A6"
    } else {
      pal[i] <- accent[((ai - 1L) %% length(accent)) + 1L]
      ai <- ai + 1L
    }
  }
  return(pal)
}

format_3 <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.3f", round(x, 3)))
}

is_integerish_column <- function(x) {
  if (!is.numeric(x)) {
    return(FALSE)
  }
  values <- stats::na.omit(x)
  if (length(values) == 0) {
    return(TRUE)
  }
  return(all(abs(values - round(values)) < sqrt(.Machine$double.eps)))
}

format_dt_table <- function(data, ...) {
  out <- DT::datatable(data, ...)
  numeric_cols <- names(data)[vapply(data, function(x) is.numeric(x) && !is_integerish_column(x), logical(1))]
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
    ggplot2::scale_color_manual(values = arm_palette(tv$Treatment)) +
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

as_display_table <- function(x) {
  if (is.data.frame(x)) {
    return(as.data.frame(x, check.names = FALSE))
  }
  if (is.atomic(x) && !is.null(names(x))) {
    return(data.frame(as.list(x), check.names = FALSE))
  }
  if (is.list(x) && !is.null(names(x))) {
    return(data.frame(as.list(x), check.names = FALSE))
  }
  return(data.frame(Value = x, check.names = FALSE))
}
