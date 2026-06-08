validation_issue <- function(type, message) {
  return(tibble::tibble(type = type, message = message))
}

validate_invivosyn_experiment <- function(
    tv, role_map, comparator_map, selected_day = NULL, min_observations = 2L) {
  errors <- list()
  warnings <- list()
  add_error <- function(message) {
    errors[[length(errors) + 1L]] <<- validation_issue("error", message)
    return(invisible(NULL))
  }
  add_warning <- function(message) {
    warnings[[length(warnings) + 1L]] <<- validation_issue("warning", message)
    return(invisible(NULL))
  }

  if (sum(role_map$role == "Vehicle") != 1L) add_error("Assign exactly one Vehicle arm.")
  if (sum(role_map$role == "Combination") < 1L) add_error("Assign at least one Combination arm.")
  if (sum(role_map$role == "Single agent") < 2L) add_error("Assign at least two Single-agent arms.")
  if (anyDuplicated(role_map$arm_id)) add_error("Each arm must have one unique role assignment.")

  combos <- role_map$arm_id[role_map$role == "Combination"]
  map_counts <- comparator_map |>
    dplyr::count(.data$combination_arm_id, name = "n")
  incomplete <- setdiff(combos, map_counts$combination_arm_id[map_counts$n >= 2L])
  if (length(incomplete) > 0) {
    add_error(paste0("Map at least two distinct comparators for: ", paste(incomplete, collapse = ", "), "."))
  }

  dup <- tv |>
    dplyr::count(.data$arm_id, .data$Mouse, .data$Day) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(dup) > 0) add_error("Treatment-Mouse-Day observations must be unique.")
  if (any(is.na(tv$Day))) add_error("Study days must be numeric and nonmissing.")
  if (any(is.na(tv$TV))) add_warning("Missing tumor-volume measurements will be excluded.")
  if (any(tv$TV < 0, na.rm = TRUE)) add_error("Tumor volumes must be nonnegative.")

  baseline_missing <- tv |>
    dplyr::summarise(
      missing = all(is.na(.data$TV0)),
      .by = c(.data$arm_id, .data$Mouse)
    ) |>
    dplyr::filter(.data$missing)
  if (nrow(baseline_missing) > 0) add_error("Every mouse must have a nonmissing baseline measurement.")
  if (any(tv$TV0 == 0, na.rm = TRUE)) add_warning("Zero baseline tumor volume makes RTV unavailable.")

  if (!is.null(selected_day) && length(combos) > 0) {
    vehicle <- role_map$arm_id[role_map$role == "Vehicle"]
    for (combo in combos) {
      comps <- comparator_map$comparator_arm_id[comparator_map$combination_arm_id == combo]
      needed <- c(vehicle, comps, combo)
      counts <- tv |>
        dplyr::filter(.data$Day == selected_day, .data$arm_id %in% needed, !is.na(.data$TV)) |>
        dplyr::count(.data$arm_id)
      missing <- setdiff(needed, counts$arm_id[counts$n >= min_observations])
      if (length(missing) > 0) {
        add_error(paste0("Day ", selected_day, " lacks enough observations for ", combo, ": ", paste(missing, collapse = ", "), "."))
      }
    }
  }

  sizes <- tv |>
    dplyr::distinct(.data$arm_id, .data$Mouse) |>
    dplyr::count(.data$arm_id)
  if (dplyr::n_distinct(sizes$n) > 1L) add_warning("Experimental arms have unequal mouse counts.")

  error_df <- dplyr::bind_rows(errors)
  warning_df <- dplyr::bind_rows(warnings)
  return(list(errors = error_df, warnings = warning_df, valid = nrow(error_df) == 0L))
}
