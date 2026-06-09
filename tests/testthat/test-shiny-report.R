test_that("report template renders from an immutable snapshot", {
  skip_if_not(rmarkdown::pandoc_available())
  app_dir <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  source(file.path(app_dir, "R", "helpers.R"), local = TRUE)
  tgi_png <- tempfile(fileext = ".png")
  grDevices::png(tgi_png, width = 400, height = 300)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
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
          synergy = c("Expected TGI" = 50, "Observed TGI" = 60, "Synergy score" = 10, "ss_lb" = 2, "ss_ub" = 18, "CI" = 1.2, "p.val.Synergy)" = 0.01, "p.val.Antagonism" = 0.99),
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
  html <- paste(readLines(rendered, warn = FALSE), collapse = "\n")
  # Package figure must be embedded as a base64 image, not printed as a raw
  # include_graphics path object.
  expect_match(html, "data:image/png;base64", fixed = TRUE)
  expect_false(grepl("knit_image_paths", html, fixed = TRUE))
})
