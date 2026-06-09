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

test_that("guided navigation wiring is present in app + helpers", {
  root <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  if (!nzchar(root)) {
    root <- testthat::test_path("..", "..", "inst", "shiny", "invivoSyn")
  }
  app <- paste(readLines(file.path(root, "app.R"), warn = FALSE), collapse = "\n")
  expect_match(app, "id = \"nav\"", fixed = TRUE)
  expect_true(all(vapply(
    c("nav_to_review", "nav_to_analyze", "nav_to_report"),
    function(x) grepl(x, app, fixed = TRUE), logical(1)
  )))
  expect_match(app, "fillable = FALSE", fixed = TRUE)
  expect_match(app, "resizable-frame", fixed = TRUE)
  expect_match(app, "split-gutter", fixed = TRUE)
  expect_match(app, ".split-pane > .bslib-card", fixed = TRUE)
  helpers <- paste(readLines(file.path(root, "R", "helpers.R"), warn = FALSE), collapse = "\n")
  expect_match(helpers, "next_button <- function", fixed = TRUE)
  expect_match(helpers, "resizable_frame <- function", fixed = TRUE)
  expect_match(helpers, "split_container <- function", fixed = TRUE)
  expect_match(helpers, "w-100", fixed = TRUE)
})
