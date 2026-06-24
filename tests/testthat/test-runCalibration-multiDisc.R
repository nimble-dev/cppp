library(testthat)
library(cppp)

## These tests exercise the generic engine directly with dummy components, so
## they need no NIMBLE model and run fast. They check that runCalibration()
## handles K discrepancies (matrix obs/sim) and stays backward compatible with
## a single discrepancy (vector obs/sim).

makeDummy <- function(discFun) {
  list(
    MCMCFun            = function(targetData, control) matrix(targetData, nrow = 5, ncol = 1,
                                                              dimnames = list(NULL, "theta")),
    simulateNewDataFun = function(thetaRow, control) as.numeric(thetaRow["theta"]),
    discFun            = discFun
  )
}

test_that("runCalibration() handles multiple discrepancies (matrix obs/sim)", {
  ## two discrepancies; sim deterministically dominates obs for disc 'a' and not 'b'
  discFun <- function(MCMCSamples, targetData, control = NULL, ...) {
    nd <- nrow(MCMCSamples)
    obs <- cbind(a = rep(0, nd), b = rep(1, nd))
    sim <- cbind(a = rep(1, nd), b = rep(0, nd))
    list(obs = obs, sim = sim, names = c("a", "b"))
  }

  s <- makeDummy(discFun)
  MCMCSamples <- matrix(1:6, ncol = 1, dimnames = list(NULL, "theta"))

  res <- runCalibration(
    MCMCSamples        = MCMCSamples,
    observedData       = 1,
    MCMCFun            = s$MCMCFun,
    simulateNewDataFun = s$simulateNewDataFun,
    discFun            = s$discFun,
    nReps              = 4
  )

  expect_named(res$CPPP, c("a", "b"))
  expect_named(res$obsPPP, c("a", "b"))
  expect_equal(dim(res$repPPP), c(4L, 2L))
  expect_equal(colnames(res$repPPP), c("a", "b"))

  ## disc 'a': sim (1) >= obs (0) always -> PPP 1; disc 'b': sim (0) >= obs (1) never -> PPP 0
  expect_equal(unname(res$obsPPP), c(1, 0))
  expect_equal(res$discrepancies$names, c("a", "b"))
})

test_that("runCalibration() stays backward compatible with a single discrepancy (vector)", {
  discFun <- function(MCMCSamples, targetData, control = NULL, ...) {
    nd <- nrow(MCMCSamples)
    list(obs = rep(0, nd), sim = rep(1, nd))
  }

  s <- makeDummy(discFun)
  MCMCSamples <- matrix(1:6, ncol = 1, dimnames = list(NULL, "theta"))

  res <- runCalibration(
    MCMCSamples        = MCMCSamples,
    observedData       = 1,
    MCMCFun            = s$MCMCFun,
    simulateNewDataFun = s$simulateNewDataFun,
    discFun            = s$discFun,
    nReps              = 4
  )

  expect_length(res$CPPP, 1L)
  expect_equal(unname(res$obsPPP), 1)       # sim always >= obs
  expect_equal(dim(res$repPPP), c(4L, 1L))
  expect_equal(colnames(res$repPPP), "discrepancy")
})

test_that("runCalibration() errors when obs/sim shapes disagree", {
  discFun <- function(MCMCSamples, targetData, control = NULL, ...) {
    nd <- nrow(MCMCSamples)
    list(obs = rep(0, nd), sim = rep(1, nd + 1))
  }

  s <- makeDummy(discFun)
  MCMCSamples <- matrix(1:6, ncol = 1, dimnames = list(NULL, "theta"))

  expect_error(
    runCalibration(
      MCMCSamples        = MCMCSamples,
      observedData       = 1,
      MCMCFun            = s$MCMCFun,
      simulateNewDataFun = s$simulateNewDataFun,
      discFun            = s$discFun,
      nReps              = 2
    ),
    "same dimensions"
  )
})
