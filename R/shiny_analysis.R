capture_package_figure <- function(plot_call, width = 1600, height = 800, res = 150) {
  file <- tempfile(pattern = "invivoSyn-figure-", fileext = ".png")
  device_open <- FALSE
  grDevices::png(filename = file, width = width, height = height, res = res)
  device_open <- TRUE
  on.exit({
    if (device_open) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  }, add = TRUE)
  value <- plot_call()
  grDevices::dev.off()
  device_open <- FALSE
  return(list(
    value = value,
    file = normalizePath(file, winslash = "/", mustWork = TRUE)
  ))
}

analysis_snapshot_id <- function(tv, role_map, comparator_map, settings) {
  return(digest::digest(list(tv, role_map, comparator_map, settings), algo = "xxhash64"))
}

analyze_combination_with_package <- function(tv, role_map, comparator_map, combo_arm_id, settings) {
  combo_tv <- build_analysis_tv(tv, role_map, comparator_map, combo_arm_id)
  combo_treatment <- lookup_treatment(role_map, combo_arm_id)

  # Only the selected metric is calculated; switching the metric re-runs analysis.
  if (identical(settings$metric, "TGI")) {
    tgi_lst <- getTGI(
      combo_tv,
      sel_day = settings$selected_day,
      tv_var = settings$tv_var,
      ci = settings$conf,
      ci_type = "bca",
      n_rep = settings$boot_n
    )
    capture <- capture_package_figure(function() {
      TGI_synergy(
        tgi_lst,
        method = settings$method,
        ci = settings$conf,
        ci_type = "bca",
        display = TRUE,
        save = FALSE
      )
    })
    efficacy <- tgi_lst$bsTGI_df
  } else {
    auc_lst <- get_mAUCr(
      combo_tv,
      sel_day = settings$end_day,
      ci = settings$conf,
      ci_type = "bca",
      nrep = settings$boot_n
    )
    capture <- capture_package_figure(function() {
      AUC_synergy(
        auc_lst,
        t = settings$auc_t,
        method = settings$method,
        boot_n = settings$boot_n,
        ci = settings$conf,
        ci_type = "bca",
        display = TRUE,
        save = FALSE,
        parallel = "snow"
      )
    })
    efficacy <- auc_lst$bsAUC_df
  }

  return(list(
    combination_arm_id = combo_arm_id,
    combination_treatment = combo_treatment,
    metric = settings$metric,
    efficacy = efficacy,
    synergy = capture$value,
    figure = capture$file
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
