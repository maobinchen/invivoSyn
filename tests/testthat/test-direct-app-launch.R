test_that("direct source launch resolves package analysis helpers", {
  app_dir <- testthat::test_path("..", "..", "inst", "shiny", "invivoSyn")
  app_env <- new.env(parent = globalenv())

  sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE)

  expect_true(exists("suggest_arm_roles", envir = app_env, inherits = FALSE))
  expect_true(exists("analyze_combinations", envir = app_env, inherits = FALSE))
  expect_identical(environment(app_env$upload_server), app_env)
  expect_identical(environment(app_env$groups_server), app_env)
})

test_that("installed app launch resolves helpers from the package namespace", {
  app_dir <- system.file("shiny", "invivoSyn", package = "invivoSyn")
  app_env <- new.env(parent = globalenv())

  expect_no_error(sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE))
  expect_true(exists("suggest_arm_roles", envir = app_env, inherits = FALSE))
  expect_true(exists("analyze_combinations", envir = app_env, inherits = FALSE))
  expect_identical(environment(app_env$upload_server), app_env)
  expect_identical(environment(app_env$groups_server), app_env)
})

test_that("direct app upload maps and parses a sample CSV", {
  app_dir <- testthat::test_path("..", "..", "inst", "shiny", "invivoSyn")
  app_env <- new.env(parent = globalenv())
  sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE)
  sample_file <- testthat::test_path("..", "..", "inst", "extdata", "test.csv")

  shiny::testServer(app_env$upload_server, {
    session$setInputs(file = list(
      name = "test.csv",
      size = file.info(sample_file)$size,
      type = "text/csv",
      datapath = sample_file
    ))
    expect_no_error(output$mapping_ui)

    session$setInputs(
      format = "auto",
      treatment = "Treatment",
      mouse = "Mouse",
      day = "",
      tv = "",
      parse = 1
    )
    expect_no_error(output$status)
  })
})

test_that("direct app upload ignores extra wide metadata columns", {
  app_dir <- testthat::test_path("..", "..", "inst", "shiny", "invivoSyn")
  app_env <- new.env(parent = globalenv())
  sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE)
  temp_csv <- tempfile(fileext = ".csv")
  dat <- data.frame(
    Treatment = c("Vehicle", "A"),
    Mouse = c("1", "2"),
    Group = c("G1", "G2"),
    `0` = c(100, 110),
    `7` = c(150, 140),
    check.names = FALSE
  )
  utils::write.csv(dat, temp_csv, row.names = FALSE)

  shiny::testServer(app_env$upload_server, {
    session$setInputs(file = list(
      name = "wide-extra.csv",
      size = file.info(temp_csv)$size,
      type = "text/csv",
      datapath = temp_csv
    ))
    session$setInputs(
      format = "wide",
      treatment = "Treatment",
      mouse = "Mouse",
      parse = 1
    )
    expect_no_error(output$status)
    expect_no_error(output$preview)
  })
})
