expected_tgi <- function(single_tgi, method = c("Bliss", "HSA")) {
  method <- match.arg(method)
  value <- switch(
    method,
    Bliss = 100 * (1 - prod(1 - single_tgi / 100)),
    HSA = max(single_tgi)
  )
  return(value)
}

expected_auc_effect <- function(single_effects, method = c("Bliss", "HSA")) {
  method <- match.arg(method)
  value <- switch(method, Bliss = prod(single_effects), HSA = min(single_effects))
  return(value)
}

bootstrap_summary <- function(scores, point_estimate, conf = 0.95) {
  alpha <- (1 - conf) / 2
  return(tibble::tibble(
    synergy_score = point_estimate,
    lb = stats::quantile(scores, alpha, na.rm = TRUE, names = FALSE),
    ub = stats::quantile(scores, 1 - alpha, na.rm = TRUE, names = FALSE),
    p_value_synergy = mean(scores <= 0, na.rm = TRUE),
    p_value_antagonism = mean(scores >= 0, na.rm = TRUE),
    interpretation = dplyr::case_when(
      stats::quantile(scores, alpha, na.rm = TRUE, names = FALSE) > 0 ~ "Synergistic",
      stats::quantile(scores, 1 - alpha, na.rm = TRUE, names = FALSE) < 0 ~ "Antagonistic",
      TRUE ~ "Inconclusive"
    )
  ))
}

resample_by_arm <- function(data) {
  return(data |>
    split(data$arm_id) |>
    purrr::map_dfr(~ dplyr::slice_sample(.x, prop = 1, replace = TRUE)))
}

analyze_tgi_combination <- function(
    tv, vehicle_arm, comparator_arms, combination_arm, selected_day,
    tv_var = c("DeltaTV", "RTV"), method = c("Bliss", "HSA"),
    boot_n = 1000L, conf = 0.95, seed = 123456L) {
  tv_var <- match.arg(tv_var)
  method <- match.arg(method)
  arms <- c(vehicle_arm, comparator_arms, combination_arm)
  day_data <- tv |>
    dplyr::filter(.data$Day == selected_day, .data$arm_id %in% arms, !is.na(.data[[tv_var]]))

  score_once <- function(x) {
    means <- x |>
      dplyr::summarise(value = mean(.data[[tv_var]]), .by = arm_id)
    vehicle <- means$value[match(vehicle_arm, means$arm_id)]
    tgi <- 100 * (1 - means$value / vehicle)
    names(tgi) <- means$arm_id
    expected <- expected_tgi(tgi[comparator_arms], method)
    observed <- tgi[[combination_arm]]
    return(c(observed = observed, expected = expected, score = observed - expected))
  }
  point <- score_once(day_data)
  set.seed(seed)
  boot <- purrr::map_dfr(seq_len(boot_n), function(iteration) {
    values <- score_once(resample_by_arm(day_data))
    return(tibble::tibble(iteration = iteration, score = values[["score"]]))
  })
  summary <- bootstrap_summary(boot$score, point[["score"]], conf) |>
    dplyr::mutate(
      combination_arm_id = combination_arm,
      metric = "TGI",
      method = method,
      observed = point[["observed"]],
      expected = point[["expected"]],
      .before = 1
    )
  return(list(summary = summary, bootstrap = boot))
}

mouse_auc <- function(data, end_day) {
  values <- data |>
    dplyr::filter(.data$Day <= end_day) |>
    dplyr::arrange(.data$Day)
  if (nrow(values) < 2L) return(NA_real_)
  span <- max(values$Day) - min(values$Day)
  if (span <= 0) return(NA_real_)
  log_tv <- log(values$TV + 1)
  auc <- pracma::trapz(values$Day, log_tv) -
    pracma::trapz(values$Day, rep(log_tv[[1]], length(log_tv)))
  return(2 * auc / (span^2))
}

analyze_auc_combination <- function(
    tv, vehicle_arm, comparator_arms, combination_arm, end_day,
    method = c("Bliss", "HSA"), boot_n = 1000L, conf = 0.95,
    seed = 123456L) {
  method <- match.arg(method)
  arms <- c(vehicle_arm, comparator_arms, combination_arm)
  filtered <- dplyr::filter(tv, .data$arm_id %in% arms)
  auc_data <- filtered |>
    split(interaction(filtered$arm_id, filtered$Mouse, drop = TRUE)) |>
    purrr::map_dfr(function(mouse_data) {
      return(tibble::tibble(
        arm_id = mouse_data$arm_id[[1]],
        Mouse = mouse_data$Mouse[[1]],
        AUC = mouse_auc(mouse_data, end_day)
      ))
    }) |>
    dplyr::filter(!is.na(.data$AUC))

  score_once <- function(x) {
    means <- x |>
      dplyr::summarise(value = mean(.data$AUC), .by = arm_id)
    vehicle <- means$value[match(vehicle_arm, means$arm_id)]
    effects <- exp((means$value - vehicle) * end_day)
    names(effects) <- means$arm_id
    expected <- expected_auc_effect(effects[comparator_arms], method)
    observed <- effects[[combination_arm]]
    return(c(observed = observed, expected = expected, score = 100 * (expected - observed)))
  }
  point <- score_once(auc_data)
  set.seed(seed)
  boot <- purrr::map_dfr(seq_len(boot_n), function(iteration) {
    values <- score_once(resample_by_arm(auc_data))
    return(tibble::tibble(iteration = iteration, score = values[["score"]]))
  })
  summary <- bootstrap_summary(boot$score, point[["score"]], conf) |>
    dplyr::mutate(
      combination_arm_id = combination_arm,
      metric = "AUC",
      method = method,
      observed = point[["observed"]],
      expected = point[["expected"]],
      .before = 1
    )
  return(list(summary = summary, bootstrap = boot))
}

rank_combination_results <- function(results) {
  return(results |>
    dplyr::arrange(.data$metric, .data$method, dplyr::desc(.data$synergy_score)) |>
    dplyr::mutate(rank = dplyr::row_number(), .by = c(metric, method)))
}

analysis_snapshot_id <- function(tv, role_map, comparator_map, settings) {
  return(digest::digest(list(tv, role_map, comparator_map, settings), algo = "xxhash64"))
}

analyze_combinations <- function(tv, role_map, comparator_map, settings) {
  vehicle <- role_map$arm_id[role_map$role == "Vehicle"][[1]]
  combos <- role_map$arm_id[role_map$role == "Combination"]
  analyses <- purrr::map(combos, function(combo) {
    comparators <- comparator_map$comparator_arm_id[
      comparator_map$combination_arm_id == combo
    ]
    if (settings$metric == "TGI") {
      return(analyze_tgi_combination(
        tv, vehicle, comparators, combo, settings$selected_day,
        settings$tv_var, settings$method, settings$boot_n, settings$conf
      ))
    }
    return(analyze_auc_combination(
      tv, vehicle, comparators, combo, settings$end_day,
      settings$method, settings$boot_n, settings$conf
    ))
  })
  summaries <- purrr::map_dfr(analyses, "summary") |>
    rank_combination_results()
  bootstraps <- purrr::map2_dfr(analyses, combos, function(x, combo) {
    return(dplyr::mutate(x$bootstrap, combination_arm_id = combo, .before = 1))
  })
  return(list(summary = summaries, bootstrap = bootstraps))
}
