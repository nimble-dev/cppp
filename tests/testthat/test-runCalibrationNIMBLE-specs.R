library(testthat)
library(nimble)
library(cppp)

## runCalibrationNIMBLE() can build the discrepancy calculator itself, from
## specifications, instead of being handed a discFun.

makeSpecTestModel <- function() {
  code <- nimbleCode({
    for (i in 1:n) y[i] ~ dnorm(mu, sd = 1)
    mu ~ dnorm(0, sd = 10)
  })

  nimbleModel(code, constants = list(n = 5),
              data = list(y = c(1, 2, 3, 4, 5)),
              inits = list(mu = 0))
}

test_that("runCalibrationNIMBLE() builds the calculator from specifications", {
  model <- makeSpecTestModel()

  ## a replicate dataset from one draw; still hand-written for now
  newDataFun <- function(thetaRow, control) {
    m  <- control$model
    nd <- m$expandNodeNames("y")
    m[["mu"]] <- thetaRow[["mu"]]
    m$simulate(nodes = nd, includeData = TRUE)
    values(m, nd)
  }

  MCMCSamples <- matrix(c(0, 0.5, 1, 1.5), ncol = 1, dimnames = list(NULL, "mu"))

  set.seed(1)
  res <- runCalibrationNIMBLE(
    model         = model,
    dataNames     = "y",
    paramNames    = "mu",
    MCMCSamples   = MCMCSamples,
    discrepancies = list(discrepancy("mean"), discrepancy("deviance")),
    simulation    = simulation("conditional"),
    simulateNewDataFun = newDataFun,
    nReps          = 2,
    MCMCcontrolRep = list(niter = 10, nburnin = 0, thin = 1),
    control        = list(model = model$newModel())
  )

  expect_named(res$CPPP, c("mean", "deviance"))
  expect_named(res$obsPPP, c("mean", "deviance"))
  expect_equal(dim(res$repPPP), c(2L, 2L))
})

test_that("runCalibrationNIMBLE() defaults to conditional simulation", {
  model <- makeSpecTestModel()

  newDataFun <- function(thetaRow, control) {
    m  <- control$model
    nd <- m$expandNodeNames("y")
    m[["mu"]] <- thetaRow[["mu"]]
    m$simulate(nodes = nd, includeData = TRUE)
    values(m, nd)
  }

  set.seed(1)
  res <- runCalibrationNIMBLE(
    model         = model,
    dataNames     = "y",
    paramNames    = "mu",
    MCMCSamples   = matrix(c(0, 1), ncol = 1, dimnames = list(NULL, "mu")),
    discrepancies = discrepancy("mean"),      # no `simulation` given
    simulateNewDataFun = newDataFun,
    nReps          = 1,
    MCMCcontrolRep = list(niter = 10, nburnin = 0, thin = 1),
    control        = list(model = model$newModel())
  )

  expect_named(res$CPPP, "mean")
})

## Nothing hand-written: the calculator and the replicate simulator both come
## from the specifications.
test_that("runCalibrationNIMBLE() runs from specifications alone", {
  model <- makeSpecTestModel()

  set.seed(1)
  res <- runCalibrationNIMBLE(
    model         = model,
    dataNames     = "y",
    paramNames    = "mu",
    MCMCSamples   = matrix(c(0, 0.5, 1), ncol = 1, dimnames = list(NULL, "mu")),
    discrepancies = list(discrepancy("mean"), discrepancy("deviance")),
    simulation    = simulation("conditional"),
    nReps          = 2,
    MCMCcontrolRep = list(niter = 10, nburnin = 0, thin = 1)
  )

  expect_named(res$CPPP, c("mean", "deviance"))
  expect_equal(dim(res$repPPP), c(2L, 2L))
})

test_that("makeSimulateNewDataFun() simulates one replicate of the right size", {
  model <- makeSpecTestModel()

  simFun <- makeSimulateNewDataFun(model$newModel(), simulation("conditional"),
                                   paramNodes = "mu")

  set.seed(1)
  rep1 <- simFun(c(mu = 100))
  rep2 <- simFun(c(mu = 100))

  expect_length(rep1, 5L)
  expect_true(all(is.finite(rep1)))
  ## drawn around mu = 100 with sd 1, and a fresh draw each time
  expect_true(all(abs(rep1 - 100) < 6))
  expect_false(identical(rep1, rep2))

  expect_error(simFun(c(sigma = 1)), "no value for")
})

## The same derived-node trap as in the calculator: with the prior on
## log(sigma), a whole-model calculate() would rebuild sigma and the replicate
## would be simulated with the wrong spread.
test_that("makeSimulateNewDataFun() uses a derived parameter as drawn", {
  code <- nimbleCode({
    for (i in 1:n) y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dflat()
    log(sigma) ~ dflat()
  })
  model <- nimbleModel(code, constants = list(n = 200), data = list(y = rep(0, 200)),
                       inits = list(mu = 0, log_sigma = 0))   # sigma = 1

  simFun <- makeSimulateNewDataFun(model$newModel(), simulation("conditional"),
                                   paramNodes = c("mu", "sigma"))

  set.seed(1)
  rep <- simFun(c(mu = 0, sigma = 10))

  ## drawn with sigma = 10, not the model's starting sigma = 1
  expect_gt(sd(rep), 5)
})

test_that("runCalibrationNIMBLE() wants one of discrepancies or discFun", {
  model <- makeSpecTestModel()
  newDataFun <- function(thetaRow, control) values(control$model, "y")

  expect_error(
    runCalibrationNIMBLE(model, dataNames = "y", paramNames = "mu",
                         MCMCSamples = matrix(0, ncol = 1, dimnames = list(NULL, "mu")),
                         simulateNewDataFun = newDataFun, nReps = 1,
                         control = list(model = model)),
    "discrepancies"
  )

  expect_error(
    runCalibrationNIMBLE(model, dataNames = "y", paramNames = "mu",
                         MCMCSamples = matrix(0, ncol = 1, dimnames = list(NULL, "mu")),
                         discrepancies = discrepancy("mean"),
                         discFun = function(...) NULL,
                         simulateNewDataFun = newDataFun, nReps = 1,
                         control = list(model = model)),
    "not both"
  )
})
