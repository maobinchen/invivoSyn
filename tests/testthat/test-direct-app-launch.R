test_that("direct source launch resolves package analysis helpers", {
  app_dir <- testthat::test_path("..", "..", "inst", "shiny", "invivoSyn")
  app_env <- new.env(parent = globalenv())

  sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE)

  expect_true(exists("suggest_arm_roles", envir = app_env, inherits = FALSE))
  expect_true(exists("analyze_combinations", envir = app_env, inherits = FALSE))
})

test_that("installed app launch resolves helpers from the package namespace", {
  app_dir <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  app_env <- new.env(parent = globalenv())

  expect_no_error(sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE))
  expect_true(exists("suggest_arm_roles", envir = app_env, inherits = FALSE))
  expect_true(exists("analyze_combinations", envir = app_env, inherits = FALSE))
})
