test_that("report template renders from an immutable snapshot", {
  skip_if_not(rmarkdown::pandoc_available())
  app_dir <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  source(file.path(app_dir, "R", "helpers.R"), local = TRUE)
  tgi_png <- tempfile(fileext = ".png")
  auc_png <- tempfile(fileext = ".png")
  grDevices::png(tgi_png, width = 100, height = 100)
  plot.new()
  grDevices::dev.off()
  grDevices::png(auc_png, width = 100, height = 100)
  plot.new()
  grDevices::dev.off()
  snapshot <- list(
    settings = list(metric = "TGI", method = "Bliss", selected_day = 7, end_day = NA, auc_t = 21),
    role_map = tibble::tibble(
      arm_id = c("Vehicle", "A", "B", "A.B"),
      Treatment = c("Vehicle", "A", "B", "A+B"),
      role = c("Vehicle", "Single agent", "Single agent", "Combination")
    ),
    comparator_map = tibble::tibble(
      combination_arm_id = "A.B",
      comparator_arm_id = c("A", "B")
    ),
    result = list(
      metric = "TGI",
      details = list(
        "A+B" = list(
          combination_arm_id = "A.B",
          combination_treatment = "A+B",
          metric = "TGI",
          efficacy = tibble::tibble(Group = c("Group 1", "Group 2"), Treatment = c("Vehicle", "A+B"), TGI = c(0, 60), std.err = c(0, 5), lb = c(0, 45), ub = c(0, 70)),
          synergy = data.frame(
            Metric = c("Bliss CI(invivoSyn)", "Bliss Synergy Score(invivoSyn)"),
            Value = c(0.18, 10), std.err = c(0.05, 4), lb = c(0.08, 2), ub = c(0.28, 18),
            p.val.Synergy = c(0.01, 0.01), p.val.Antagonism = c(0.99, 0.99),
            check.names = FALSE, stringsAsFactors = FALSE
          ),
          figure = normalizePath(tgi_png, winslash = "/", mustWork = TRUE)
        )
      )
    ),
    tv = tibble::tibble(
      arm_id = rep(c("Vehicle", "A.B"), each = 2),
      Treatment = rep(c("Vehicle", "A+B"), each = 2),
      Mouse = rep("1", 4), Day = rep(c(0, 7), 2), TV = c(100, 200, 100, 120)
    )
  )
  source_filename <- "synthetic.csv"
  output <- tempfile(fileext = ".html")
  rendered <- rmarkdown::render(
    file.path(app_dir, "report_template.Rmd"),
    output_file = basename(output), output_dir = dirname(output),
    envir = environment(), quiet = TRUE
  )
  expect_true(file.exists(rendered))
})
