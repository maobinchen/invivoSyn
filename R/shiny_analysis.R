package_tgi_bootstrap <- function(tgi_lst, method = c("Bliss", "HSA")) {
  method <- match.arg(method)
  roles <- tgi_lst$roles
  boot_matrix <- as.data.frame(tgi_lst$bsTGI_r$t)
  colnames(boot_matrix) <- c(as.character(roles$single_grps), "Combo")
  single_matrix <- as.matrix(boot_matrix[, seq_along(roles$single_grps), drop = FALSE])
  expected <- switch(
    method,
    Bliss = 100 * (1 - apply(1 - single_matrix / 100, 1, prod)),
    HSA = apply(single_matrix, 1, max)
  )
  observed <- 100 * boot_matrix$Combo
  return(tibble::tibble(
    iteration = seq_len(nrow(boot_matrix)),
    observed = observed,
    expected = expected,
    score = observed - expected
  ))
}

package_auc_bootstrap <- function(
    auc_lst, end_day, method = c("Bliss", "HSA"), boot_n = 1000L,
    seed = 123456L) {
  method <- match.arg(method)
  roles <- auc_lst$roles
  auc_mouse <- as.data.frame(auc_lst$auc_mouse)
  set.seed(seed)
  boot_idx <- replicate(
    boot_n,
    sample.int(nrow(auc_mouse), size = nrow(auc_mouse), replace = TRUE),
    simplify = FALSE
  )
  boot_scores <- purrr::map_dfr(seq_along(boot_idx), function(iteration) {
    values <- getFromNamespace("bs_AUC_synergy", "invivoSyn")(
      auc_mouse = auc_mouse,
      t = end_day,
      method = method,
      roles = roles,
      idx = boot_idx[[iteration]]
    )
    ci_value <- unname(values[["CI"]])
    synergy_score <- unname(values[["Synergy_score"]])
    if (isTRUE(all.equal(ci_value, 1))) {
      expected <- NA_real_
      observed <- NA_real_
    } else {
      expected <- (synergy_score / 100) / (1 - ci_value)
      observed <- ci_value * expected
    }
    return(tibble::tibble(
      iteration = iteration,
      observed = observed,
      expected = expected,
      score = synergy_score
    ))
  })
  return(boot_scores)
}

interpret_score <- function(lb, ub) {
  if (is.na(lb) || is.na(ub)) {
    return(NA_character_)
  }
  if (lb > 0) {
    return("Synergistic")
  }
  if (ub < 0) {
    return("Antagonistic")
  }
  return("Inconclusive")
}

summarize_tgi_result <- function(combo_arm_id, combo_treatment, method, tgi_syn, boot_scores) {
  p_syn_name <- grep("^p\\.val\\.Synergy", names(tgi_syn), value = TRUE)[1]
  summary <- tibble::tibble(
    combination_arm_id = combo_arm_id,
    combination_treatment = combo_treatment,
    metric = "TGI",
    method = method,
    observed = unname(tgi_syn[["Observed TGI"]]),
    expected = unname(tgi_syn[["Expected TGI"]]),
    synergy_score = unname(tgi_syn[["Synergy score"]]),
    lb = unname(tgi_syn[["ss_lb"]]),
    ub = unname(tgi_syn[["ss_ub"]]),
    ci_value = unname(tgi_syn[["CI"]]),
    p_value_synergy = unname(tgi_syn[[p_syn_name]]),
    p_value_antagonism = unname(tgi_syn[["p.val.Antagonism"]])
  ) |>
    dplyr::mutate(
      interpretation = interpret_score(.data$lb, .data$ub)
    )
  bootstrap <- dplyr::mutate(
    boot_scores,
    combination_arm_id = combo_arm_id,
    combination_treatment = combo_treatment,
    metric = "TGI",
    .before = 1
  )
  return(list(summary = summary, bootstrap = bootstrap))
}

summarize_auc_result <- function(combo_arm_id, combo_treatment, method, auc_syn, boot_scores) {
  synergy_row <- auc_syn[auc_syn$Metric == "Synergy_score", , drop = FALSE]
  ci_row <- auc_syn[auc_syn$Metric == "CI", , drop = FALSE]
  summary <- tibble::tibble(
    combination_arm_id = combo_arm_id,
    combination_treatment = combo_treatment,
    metric = "AUC",
    method = method,
    observed = mean(boot_scores$observed, na.rm = TRUE),
    expected = mean(boot_scores$expected, na.rm = TRUE),
    synergy_score = synergy_row$Value[[1]],
    lb = synergy_row$lb[[1]],
    ub = synergy_row$ub[[1]],
    ci_value = ci_row$Value[[1]],
    p_value_synergy = synergy_row$`p.val.Synergy`[[1]],
    p_value_antagonism = synergy_row$`p.val.Antagonism`[[1]]
  ) |>
    dplyr::mutate(
      interpretation = interpret_score(.data$lb, .data$ub)
    )
  bootstrap <- dplyr::mutate(
    boot_scores,
    combination_arm_id = combo_arm_id,
    combination_treatment = combo_treatment,
    metric = "AUC",
    .before = 1
  )
  return(list(summary = summary, bootstrap = bootstrap))
}

rank_combination_results <- function(results) {
  return(results |>
    dplyr::arrange(.data$metric, .data$method, dplyr::desc(.data$synergy_score)) |>
    dplyr::mutate(rank = dplyr::row_number(), .by = c(metric, method)))
}

analysis_snapshot_id <- function(tv, role_map, comparator_map, settings) {
  return(digest::digest(list(tv, role_map, comparator_map, settings), algo = "xxhash64"))
}

analyze_combination_with_package <- function(tv, role_map, comparator_map, combo_arm_id, settings) {
  combo_tv <- build_analysis_tv(tv, role_map, comparator_map, combo_arm_id)
  combo_treatment <- lookup_treatment(role_map, combo_arm_id)
  if (settings$metric == "TGI") {
    tgi_lst <- getTGI(
      combo_tv,
      sel_day = settings$selected_day,
      tv_var = settings$tv_var,
      ci = settings$conf,
      n_rep = settings$boot_n
    )
    tgi_syn <- TGI_synergy(
      tgi_lst,
      method = settings$method,
      ci = settings$conf,
      display = FALSE,
      save = FALSE
    )
    boot_scores <- package_tgi_bootstrap(tgi_lst, settings$method)
    return(summarize_tgi_result(combo_arm_id, combo_treatment, settings$method, tgi_syn, boot_scores))
  }
  auc_lst <- get_mAUCr(
    combo_tv,
    sel_day = settings$end_day,
    ci = settings$conf,
    nrep = settings$boot_n
  )
  auc_syn <- AUC_synergy(
    auc_lst,
    t = settings$end_day,
    method = settings$method,
    boot_n = settings$boot_n,
    ci = settings$conf,
    save = FALSE,
    display = FALSE
  )
  boot_scores <- package_auc_bootstrap(
    auc_lst,
    end_day = settings$end_day,
    method = settings$method,
    boot_n = settings$boot_n
  )
  return(summarize_auc_result(combo_arm_id, combo_treatment, settings$method, auc_syn, boot_scores))
}

analyze_combinations <- function(tv, role_map, comparator_map, settings) {
  combos <- role_map$arm_id[role_map$role == "Combination"]
  analyses <- purrr::map(
    combos,
    ~ analyze_combination_with_package(tv, role_map, comparator_map, .x, settings)
  )
  summaries <- purrr::map_dfr(analyses, "summary") |>
    rank_combination_results()
  bootstraps <- purrr::map_dfr(analyses, "bootstrap")
  return(list(summary = summaries, bootstrap = bootstraps))
}
