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

test_that("expected effects support higher-order combinations", {
  expect_equal(invivoSyn:::expected_tgi(c(20, 30, 40), "HSA"), 40)
  expect_equal(invivoSyn:::expected_tgi(c(20, 30), "Bliss"), 44)
  expect_equal(invivoSyn:::expected_auc_effect(c(0.8, 0.7), "Bliss"), 0.56)
})

test_that("comparator suggestions use combination names", {
  roles <- invivoSyn:::suggest_arm_roles(c("Vehicle", "DrugA", "DrugB", "DrugA+DrugB"))
  mapping <- invivoSyn:::suggest_comparator_map(roles)
  expect_equal(nrow(mapping), 2)
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

test_that("multiple combinations are analyzed independently and ranked", {
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
    tv_var = "RTV", conf = 0.95, boot_n = 25L
  )
  result <- invivoSyn:::analyze_combinations(tv, roles, mapping, settings)
  expect_equal(nrow(result$summary), 2)
  expect_equal(sort(result$summary$rank), 1:2)
  expect_equal(sort(unique(result$bootstrap$combination_arm_id)), sort(combos))

  settings$metric <- "AUC"
  auc_result <- invivoSyn:::analyze_combinations(tv, roles, mapping, settings)
  expect_equal(nrow(auc_result$summary), 2)
  expect_equal(sort(auc_result$summary$rank), 1:2)
})
