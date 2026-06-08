make_arm_id <- function(x) {
  ids <- make.unique(make.names(as.character(x)))
  return(ids)
}

lookup_treatment <- function(role_map, arm_id) {
  value <- role_map$Treatment[match(arm_id, role_map$arm_id)]
  return(as.character(value[[1]]))
}

lookup_arm_id <- function(role_map, treatment) {
  value <- role_map$arm_id[match(treatment, role_map$Treatment)]
  return(as.character(value[[1]]))
}

suggest_tv_columns <- function(data) {
  nm <- names(data)
  lower <- tolower(nm)
  pick <- function(patterns) {
    idx <- which(vapply(
      lower,
      function(x) any(vapply(patterns, grepl, logical(1), x = x)),
      logical(1)
    ))
    if (length(idx) == 0) return(NA_character_)
    return(nm[[idx[[1]]]])
  }
  return(list(
    treatment = pick(c("^treatment$", "^group$", "arm")),
    mouse = pick(c("^mouse$", "subject", "animal")),
    day = pick(c("^day$", "time")),
    tv = pick(c("^tv$", "tumou?r.?volume", "volume"))
  ))
}

suggest_arm_roles <- function(treatments) {
  treatment <- unique(as.character(treatments))
  low <- tolower(treatment)
  role <- dplyr::case_when(
    grepl("vehicle|control|placebo", low) ~ "Vehicle",
    grepl("\\+|combo|combination", low) ~ "Combination",
    TRUE ~ "Single agent"
  )
  return(tibble::tibble(
    arm_id = make_arm_id(treatment),
    Treatment = treatment,
    role = role
  ))
}

extract_arm_tokens <- function(label) {
  cleaned <- tolower(label)
  cleaned <- gsub("^group\\s*[0-9]+\\s*[,._ -]*", "", cleaned)
  cleaned <- gsub("[()]", " ", cleaned)
  cleaned <- gsub("[^a-z0-9]+", " ", cleaned)
  tokens <- unlist(strsplit(cleaned, "\\s+"))
  tokens <- tokens[nzchar(tokens)]
  return(tokens)
}

build_comparator_map <- function(role_map, comparator_map) {
  combinations <- role_map |>
    dplyr::filter(.data$role == "Combination") |>
    dplyr::pull(.data$arm_id)
  singles <- role_map |>
    dplyr::filter(.data$role == "Single agent") |>
    dplyr::pull(.data$arm_id)

  out <- comparator_map |>
    dplyr::transmute(
      combination_arm_id = as.character(.data$combination_arm_id),
      comparator_arm_id = as.character(.data$comparator_arm_id)
    ) |>
    dplyr::distinct()

  invalid_combo <- setdiff(out$combination_arm_id, combinations)
  invalid_single <- setdiff(out$comparator_arm_id, singles)
  if (length(invalid_combo) > 0 || length(invalid_single) > 0) {
    rlang::abort("Comparator mappings must connect Combination arms to Single-agent arms.")
  }
  return(out)
}

build_combo_role_args <- function(role_map, comparator_map, combo_arm_id) {
  vehicle_arm_id <- role_map$arm_id[role_map$role == "Vehicle"]
  comparator_arm_ids <- comparator_map$comparator_arm_id[
    comparator_map$combination_arm_id == combo_arm_id
  ]
  return(list(
    vehicle = lookup_treatment(role_map, vehicle_arm_id[[1]]),
    singles = unname(vapply(comparator_arm_ids, lookup_treatment, character(1), role_map = role_map)),
    combo = lookup_treatment(role_map, combo_arm_id)
  ))
}

suggest_comparator_map <- function(role_map) {
  singles <- role_map |>
    dplyr::filter(.data$role == "Single agent")
  combinations <- role_map |>
    dplyr::filter(.data$role == "Combination")
  suggested <- purrr::map_dfr(seq_len(nrow(combinations)), function(i) {
    combo <- combinations[i, ]
    combo_tokens <- extract_arm_tokens(combo$Treatment)
    matched <- singles$arm_id[vapply(
      singles$Treatment,
      function(single_name) {
        single_tokens <- extract_arm_tokens(single_name)
        any(single_tokens %in% combo_tokens)
      },
      logical(1)
    )]
    return(tibble::tibble(
      combination_arm_id = combo$arm_id,
      comparator_arm_id = matched
    ))
  })
  return(suggested)
}

latest_common_day <- function(tv, arm_ids, min_observations = 2L) {
  counts <- tv |>
    dplyr::filter(.data$arm_id %in% arm_ids, !is.na(.data$TV)) |>
    dplyr::count(.data$arm_id, .data$Day, name = "n")
  valid_days <- counts |>
    dplyr::filter(.data$n >= min_observations) |>
    dplyr::summarise(n_arms = dplyr::n_distinct(.data$arm_id), .by = Day) |>
    dplyr::filter(.data$n_arms == length(unique(arm_ids))) |>
    dplyr::arrange(.data$Day)
  if (nrow(valid_days) == 0) {
    return(NA_real_)
  }
  return(max(valid_days$Day))
}

build_analysis_tv <- function(tv, role_map, comparator_map, combo_arm_id) {
  role_args <- build_combo_role_args(role_map, comparator_map, combo_arm_id)
  out <- set_roles(
    tv,
    vehicle = role_args$vehicle,
    singles = role_args$singles,
    combo = role_args$combo
  )
  roles_full <- attr(out, "roles")
  out <- dplyr::mutate(
    out,
    Mouse = paste(as.character(.data$Treatment), .data$Mouse, sep = "::")
  )
  attr(out, "roles") <- roles_full
  attr(out, "analysis_combo_arm_id") <- combo_arm_id
  return(out)
}

normalize_tv_long <- function(data, treatment_col, mouse_col, day_col, tv_col) {
  required <- c(treatment_col, mouse_col, day_col, tv_col)
  if (any(!required %in% names(data))) {
    rlang::abort("One or more mapped columns are absent from the uploaded data.")
  }
  out <- data |>
    dplyr::transmute(
      Treatment = as.character(.data[[treatment_col]]),
      Mouse = as.character(.data[[mouse_col]]),
      Day = suppressWarnings(as.numeric(.data[[day_col]])),
      TV = suppressWarnings(as.numeric(.data[[tv_col]]))
    )
  arms <- tibble::tibble(
    Treatment = unique(out$Treatment),
    arm_id = make_arm_id(unique(out$Treatment)),
    Group = paste0("Group ", seq_along(unique(out$Treatment)))
  )
  out <- out |>
    dplyr::left_join(arms, by = "Treatment") |>
    dplyr::mutate(
      baseline_day = min(.data$Day, na.rm = TRUE),
      TV0 = .data$TV[match(.data$baseline_day[[1]], .data$Day)],
      .by = c(arm_id, Mouse)
    ) |>
    dplyr::mutate(
      RTV = dplyr::if_else(.data$TV0 == 0, NA_real_, .data$TV / .data$TV0),
      DeltaTV = .data$TV - .data$TV0,
      logTV = log(.data$TV + 1)
    ) |>
    dplyr::select(
      "arm_id", "Group", "Treatment", "Mouse", "Day",
      "TV", "TV0", "RTV", "DeltaTV", "logTV"
    )
  return(out)
}

normalize_tv_wide <- function(data, treatment_col, mouse_col) {
  candidate_cols <- setdiff(names(data), c(treatment_col, mouse_col))
  day_cols <- candidate_cols[!is.na(suppressWarnings(as.numeric(candidate_cols)))]
  if (length(day_cols) == 0L) {
    rlang::abort("No numeric day columns were found after excluding treatment and mouse columns.")
  }
  long <- data |>
    tidyr::pivot_longer(
      dplyr::all_of(day_cols),
      names_to = "Day",
      values_to = "TV"
    )
  return(normalize_tv_long(
    long,
    treatment_col = treatment_col,
    mouse_col = mouse_col,
    day_col = "Day",
    tv_col = "TV"
  ))
}
