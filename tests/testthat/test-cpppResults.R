# tests/testthat/test-cpppResult.R
test_that("newCpppResult creates a cpppResult with expected fields", {
  x <- newCpppResult(CPPP = 0.4, repPPP = runif(5), obsPPP = 0.3)
  expect_s3_class(x, "cpppResult")
  expect_true(is.list(x))
  expect_equal(length(x$repPPP), 5)
  expect_true(x$CPPP >= 0 && x$CPPP <= 1)
})

test_that("newCpppResult enforces bounds when finite", {
  expect_error(newCpppResult(CPPP = 1.2), "in \\[0,1\\]")
  expect_error(newCpppResult(repPPP = c(-0.1, 0.2)), "in \\[0,1\\]")
  expect_error(newCpppResult(obsPPP = c(0.5, 1.1)), "in \\[0,1\\]")
})
