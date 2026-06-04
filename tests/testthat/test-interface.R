# tests/testthat/test-interface.R


test_that("transferAutocorrelation exists when implemented", {
  expect_true(
    exists("transferAutocorrelation", envir = asNamespace("cppp"), inherits = FALSE),
    info = "transferAutocorrelation() is not yet defined or exported"
  )
})

test_that("runCalibration exists", {
  expect_true(
    exists("runCalibration", envir = asNamespace("cppp"), inherits = FALSE),
    info = "runCalibration() not found in package namespace"
  )
})
