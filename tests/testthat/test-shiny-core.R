test_that("normalization and role suggestions support multiple combinations", {
  dat <- data.frame(
    trt = rep(c("Vehicle", "A", "B", "A+B low", "A+B high"), each = 4),
    mouse = rep(rep(1:2, each = 2), 5),
    day = rep(c(0, 7), 10),
    tv = seq_len(20)
  )
  tv <- invivoSyn:::normalize_tv_long(dat, "trt", "mouse", "day", "tv")
  roles <- invivoSyn:::suggest_arm_roles(unique(dat$trt))
  expect_true(all(c("arm_id", "TV0", "RTV", "DeltaTV", "logTV") %in% names(tv)))
  expect_equal(sum(roles$role == "Combination"), 2)
})

test_that("wide normalization ignores extra non-day metadata columns", {
  dat <- data.frame(
    Treatment = c("Vehicle", "A"),
    Mouse = c("1", "2"),
    Group = c("G1", "G2"),
    `0` = c(100, 110),
    `7` = c(150, 140),
    check.names = FALSE
  )
  tv <- invivoSyn:::normalize_tv_wide(dat, "Treatment", "Mouse")
  expect_equal(sort(unique(tv$Day)), c(0, 7))
  expect_false(any(grepl("Group\\.", names(tv), fixed = FALSE)))
})

test_that("comparator suggestions ignore group prefixes", {
  roles <- invivoSyn:::suggest_arm_roles(c(
    "Group 01,Vehicle",
    "Group 04,TNO155",
    "Group 05,Trametinib",
    "Group 06,TNO155+Trametinib"
  ))
  mapping <- invivoSyn:::suggest_comparator_map(roles)
  expect_equal(sort(mapping$comparator_arm_id), sort(c("Group.04.TNO155", "Group.05.Trametinib")))
})

test_that("build_analysis_tv re-tags package roles per combination", {
  tv <- invivoSyn:::normalize_tv_long(
    data.frame(
      trt = rep(c("Vehicle", "A", "B", "A+B low", "A+B high"), each = 4),
      mouse = rep(rep(1:2, each = 2), 5),
      day = rep(c(0, 7), 10),
      tv = seq_len(20)
    ),
    "trt", "mouse", "day", "tv"
  )
  roles <- invivoSyn:::suggest_arm_roles(unique(tv$Treatment))
  singles <- roles$arm_id[roles$role == "Single agent"]
  combos <- roles$arm_id[roles$role == "Combination"]
  mapping <- tidyr::crossing(combination_arm_id = combos, comparator_arm_id = singles)

  combo_tv <- invivoSyn:::build_analysis_tv(tv, roles, mapping, combos[[1]])
  combo_roles <- invivoSyn::get_roles(combo_tv)
  expect_equal(combo_roles$vehicle, "Vehicle")
  expect_equal(sort(combo_roles$singles), sort(c("A", "B")))
  expect_equal(combo_roles$combo, "A+B low")
})

test_that("validation accepts exact mappings for multiple combinations", {
  tv <- invivoSyn:::normalize_tv_long(
    data.frame(
      trt = rep(c("Vehicle", "A", "B", "A+B 1", "A+B 2"), each = 4),
      mouse = rep(rep(1:2, each = 2), 5),
      day = rep(c(0, 7), 10),
      tv = rep(c(100, 200, 100, 150, 100, 160, 100, 120, 100, 110), each = 2)
    ),
    "trt", "mouse", "day", "tv"
  )
  roles <- invivoSyn:::suggest_arm_roles(unique(tv$Treatment))
  singles <- roles$arm_id[roles$role == "Single agent"]
  combos <- roles$arm_id[roles$role == "Combination"]
  mapping <- tidyr::crossing(combination_arm_id = combos, comparator_arm_id = singles)
  result <- invivoSyn:::validate_invivosyn_experiment(tv, roles, mapping, selected_day = 7)
  expect_true(result$valid)
})

test_that("multiple combinations are analyzed with package-backed TGI orchestration", {
  treatments <- c("Vehicle", "A", "B", "A+B low", "A+B high")
  data <- expand.grid(
    trt = treatments, mouse = as.character(1:4), day = c(0, 7, 14),
    stringsAsFactors = FALSE
  )
  effects <- c(Vehicle = 1, A = 0.75, B = 0.7, "A+B low" = 0.45, "A+B high" = 0.3)
  data$tv <- 100 * exp(0.07 * data$day) * effects[data$trt] + as.numeric(data$mouse)
  tv <- invivoSyn:::normalize_tv_long(data, "trt", "mouse", "day", "tv")
  roles <- invivoSyn:::suggest_arm_roles(treatments)
  singles <- roles$arm_id[roles$role == "Single agent"]
  combos <- roles$arm_id[roles$role == "Combination"]
  mapping <- tidyr::crossing(combination_arm_id = combos, comparator_arm_id = singles)
  settings <- list(
    metric = "TGI", method = "Bliss", selected_day = 14, end_day = 14,
    tv_var = "RTV", conf = 0.95, boot_n = 20L
  )
  result <- suppressWarnings(invivoSyn:::analyze_combinations(tv, roles, mapping, settings))
  expect_equal(nrow(result$summary), 2)
  expect_equal(sort(result$summary$rank), 1:2)
  expect_true(all(c(
    "combination_arm_id", "combination_treatment", "metric", "method",
    "observed", "expected", "synergy_score", "lb", "ub",
    "p_value_synergy", "p_value_antagonism", "interpretation"
  ) %in% names(result$summary)))
  expect_true(all(c("combination_treatment", "score") %in% names(result$bootstrap)))
})

test_that("multiple combinations are analyzed with package-backed AUC orchestration", {
  treatments <- c("Vehicle", "A", "B", "A+B low", "A+B high")
  data <- expand.grid(
    trt = treatments, mouse = as.character(1:4), day = c(0, 7, 14),
    stringsAsFactors = FALSE
  )
  effects <- c(Vehicle = 1, A = 0.75, B = 0.7, "A+B low" = 0.45, "A+B high" = 0.3)
  data$tv <- 100 * exp(0.07 * data$day) * effects[data$trt] + as.numeric(data$mouse)
  tv <- invivoSyn:::normalize_tv_long(data, "trt", "mouse", "day", "tv")
  roles <- invivoSyn:::suggest_arm_roles(treatments)
  singles <- roles$arm_id[roles$role == "Single agent"]
  combos <- roles$arm_id[roles$role == "Combination"]
  mapping <- tidyr::crossing(combination_arm_id = combos, comparator_arm_id = singles)
  settings <- list(
    metric = "AUC", method = "HSA", selected_day = 14, end_day = 14,
    tv_var = "RTV", conf = 0.95, boot_n = 20L
  )
  result <- suppressWarnings(invivoSyn:::analyze_combinations(tv, roles, mapping, settings))
  expect_equal(nrow(result$summary), 2)
  expect_true(all(result$summary$metric == "AUC"))
  expect_true(all(result$summary$method == "HSA"))
})

test_that("latest_common_day finds the last analyzable day", {
  tv <- tibble::tibble(
    arm_id = c("Vehicle", "Vehicle", "A", "A", "Combo", "Combo", "Combo"),
    Mouse = c("1", "1", "1", "1", "1", "2", "1"),
    Day = c(0, 14, 0, 14, 0, 0, 7),
    TV = c(100, 150, 100, 120, 100, 110, 115)
  )
  expect_equal(invivoSyn:::latest_common_day(tv, c("Vehicle", "A", "Combo"), min_observations = 1L), 0)
})

test_that("summary endpoint table is wide and uses mean+/-sd(n) format", {
  source(testthat::test_path("..", "..", "inst", "shiny", "invivoSyn", "R", "helpers.R"), local = TRUE)
  tv <- tibble::tibble(
    arm_id = c("A", "A", "A", "A"),
    Group = c("Group 1", "Group 1", "Group 1", "Group 1"),
    Treatment = c("Vehicle", "Vehicle", "Vehicle", "Vehicle"),
    Mouse = c("1", "2", "1", "2"),
    Day = c(0, 0, 7, 7),
    TV = c(100, 110, 150, 170),
    TV0 = c(100, 110, 100, 110),
    RTV = c(1, 1, 1.5, 1.545),
    DeltaTV = c(0, 0, 50, 60),
    logTV = log(c(101, 111, 151, 171))
  )
  summary <- summarize_endpoint_wide(tv, "TV")
  expect_true(all(c("Group", "Treatment", "0", "7") %in% names(summary)))
  expect_match(summary$`0`[[1]], "105\\.000\\+/-7\\.071\\(2\\)")
  expect_match(summary$`7`[[1]], "160\\.000\\+/-14\\.142\\(2\\)")
})

test_that("AUC bootstrap extractor supports unnamed vectors", {
  extracted <- invivoSyn:::extract_auc_bootstrap_values(c(0.8, 12.5))
  expect_equal(extracted$ci_value, 0.8)
  expect_equal(extracted$synergy_score, 12.5)
})
