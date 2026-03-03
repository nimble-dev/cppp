# tests/testthat/test-interface.R


test_that("transferAutocorrelation exists when implemented", {
  expect_true(
    "transferAutocorrelation" %in% ls("package:cpppPrototype"),
    info = "transferAutocorrelation() is not yet defined or exported"
  )
})

test_that("runCalibration exists", {
  expect_true(
    "runCalibration" %in% ls("package:cpppPrototype"),
    info = "runCalibration() not found in package namespace"
  )
})
