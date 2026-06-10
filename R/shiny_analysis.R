analysis_snapshot_id <- function(tv, role_map, comparator_map, settings) {
  return(digest::digest(list(tv, role_map, comparator_map, settings), algo = "xxhash64"))
}

analyze_combination_with_package <- function(tv, role_map, comparator_map, combo_arm_id, settings) {
  combo_tv <- build_analysis_tv(tv, role_map, comparator_map, combo_arm_id)
  combo_treatment <- lookup_treatment(role_map, combo_arm_id)

  # Only the selected metric is calculated; switching the metric re-runs analysis.
  # The synergy functions write their own 300-dpi figure via save=TRUE/file=; use
  # that file directly rather than re-capturing the printed plot.
  fig_base <- tempfile(pattern = "invivoSyn-figure-")
  auc_lst <- get_mAUCr(
    combo_tv,
    sel_day = settings$end_day,
    ci = settings$conf,
    ci_type = "bca",
    nrep = settings$boot_n
  )
  value <- AUC_synergy(
    auc_lst,
    t = settings$auc_t,
    method = settings$method,
    boot_n = settings$boot_n,
    ci = settings$conf,
    ci_type = "bca",
    display = FALSE,
    save = TRUE,
    file = fig_base,
    parallel = "snow"
  )
  efficacy <- auc_lst$bsAUC_df

  return(list(
    combination_arm_id = combo_arm_id,
    combination_treatment = combo_treatment,
    metric = settings$metric,
    efficacy = efficacy,
    synergy = value,
    figure = normalizePath(paste0(fig_base, ".png"), winslash = "/", mustWork = TRUE)
  ))
}

analyze_combinations <- function(tv, role_map, comparator_map, settings) {
  combos <- role_map$arm_id[role_map$role == "Combination"]
  analyses <- purrr::map(
    combos,
    ~ analyze_combination_with_package(tv, role_map, comparator_map, .x, settings)
  )
  details <- purrr::set_names(analyses, purrr::map_chr(analyses, "combination_treatment"))
  return(list(
    metric = settings$metric,
    details = details
  ))
}
