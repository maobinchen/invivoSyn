test_that("bundled Shiny app has required entry points", {
  root <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  if (!nzchar(root)) {
    root <- testthat::test_path("..", "..", "inst", "shiny", "invivoSyn")
  }
  expect_true(file.exists(file.path(root, "app.R")))
  expect_true(file.exists(file.path(root, "report_template.Rmd")))
  expect_true(all(file.exists(file.path(root, "R", c(
    "helpers.R", "mod_upload.R", "mod_groups.R", "mod_review.R",
    "mod_analysis.R", "mod_report.R"
  )))))
})
