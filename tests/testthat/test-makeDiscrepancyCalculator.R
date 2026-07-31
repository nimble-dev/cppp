library(testthat)
library(nimble)
library(cppp)

## A small model reused by the tests below.
makeTestModel <- function() {
  code <- nimbleCode({
    for (i in 1:n) {
      y[i] ~ dnorm(mu, sd = 1)
    }
    mu ~ dnorm(0, sd = 10)
  })

  nimbleModel(
    code = code,
    constants = list(n = 4),
    data = list(y = c(1, 2, 3, 4)),
    inits = list(mu = 0)
  )
}

test_that("makeDiscrepancyCalculator() gives one column per discrepancy", {
  model <- makeTestModel()

  discFun <- makeDiscrepancyCalculator(
    model,
    discrepancies = list(discrepancy("mean"), discrepancy("deviance")),
    simulation    = simulation("conditional"),
    paramNodes    = "mu"
  )

  MCMCSamples <- matrix(c(0, 0.5, 1), ncol = 1, dimnames = list(NULL, "mu"))
  out <- discFun(MCMCSamples, targetData = c(1, 2, 3, 4))

  expect_equal(dim(out$obs), c(3L, 2L))
  expect_equal(dim(out$sim), c(3L, 2L))
  expect_equal(colnames(out$obs), c("mean", "deviance"))
  expect_equal(colnames(out$sim), c("mean", "deviance"))
  expect_true(all(is.finite(out$sim)))
})

## The observed side is a plain function of the dataset and the parameters, so
## we can say what the numbers must be.
test_that("makeDiscrepancyCalculator() computes the observed side correctly", {
  model <- makeTestModel()
  y <- c(1, 2, 3, 4)

  discFun <- makeDiscrepancyCalculator(
    model,
    discrepancies = list(discrepancy("mean"), discrepancy("deviance")),
    simulation    = simulation("conditional"),
    paramNodes    = "mu"
  )

  mu <- c(0, 0.5, 1)
  MCMCSamples <- matrix(mu, ncol = 1, dimnames = list(NULL, "mu"))
  out <- discFun(MCMCSamples, targetData = y)

  ## mean(y) does not involve mu at all
  expect_equal(unname(out$obs[, "mean"]), rep(mean(y), 3))

  ## deviance = -2 * log-likelihood of y under each mu
  expected <- vapply(mu, function(m) -2 * sum(dnorm(y, m, 1, log = TRUE)), numeric(1))
  expect_equal(unname(out$obs[, "deviance"]), expected)
})

## A parameter can be a derived node: with log(sigma) ~ dflat(), `sigma` is
## computed from `log_sigma`. Recalculating the whole model would rebuild
## `sigma` and throw away the value the draw gave us, leaving every draw with
## the same sigma.
test_that("makeDiscrepancyCalculator() keeps a derived parameter at its drawn value", {
  code <- nimbleCode({
    for (i in 1:n) y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dflat()
    log(sigma) ~ dflat()
  })

  model <- nimbleModel(code, constants = list(n = 3), data = list(y = c(1, 2, 3)),
                       inits = list(mu = 0, log_sigma = 2))
  expect_true("sigma" %in% model$getNodeNames(determOnly = TRUE))

  y <- c(1, 2, 3)
  discFun <- makeDiscrepancyCalculator(model, discrepancy("deviance"),
                                       simulation("conditional"),
                                       paramNodes = c("mu", "sigma"))

  MCMCSamples <- cbind(mu = c(0, 0), sigma = c(1, 4))
  out <- discFun(MCMCSamples, targetData = y)

  ## deviance depends on sigma, so a frozen sigma would give two equal rows
  expected <- c(-2 * sum(dnorm(y, 0, 1, log = TRUE)),
                -2 * sum(dnorm(y, 0, 4, log = TRUE)))
  expect_equal(unname(out$obs[, "deviance"]), expected)
})

test_that("makeDiscrepancyCalculator() puts the model back as it found it", {
  model <- makeTestModel()

  discFun <- makeDiscrepancyCalculator(model, discrepancy("mean"),
                             simulation("conditional"), paramNodes = "mu")

  yBefore  <- values(model, "y")
  muBefore <- values(model, "mu")

  discFun(matrix(c(0, 1), ncol = 1, dimnames = list(NULL, "mu")),
          targetData = c(9, 9, 9, 9))

  expect_equal(values(model, "y"), yBefore)
  expect_equal(values(model, "mu"), muBefore)
})

test_that("makeDiscrepancyCalculator() complains instead of returning wrong numbers", {
  model <- makeTestModel()

  ## a discrepancy reading data the simulation never writes
  expect_error(
    makeDiscrepancyCalculator(model,
                    discrepancy("mean"),                       # all of y
                    simulation("conditional", dataNodes = "y[1:2]"),
                    paramNodes = "mu"),
    "does not set"
  )

  ## two columns that would both be called "mean"
  expect_error(
    makeDiscrepancyCalculator(model, list(discrepancy("mean"), discrepancy("mean")),
                    simulation("conditional"), paramNodes = "mu"),
    "own name"
  )

  discFun <- makeDiscrepancyCalculator(model, discrepancy("mean"),
                             simulation("conditional"), paramNodes = "mu")

  ## draws with no column for mu
  expect_error(
    discFun(matrix(0, ncol = 1, dimnames = list(NULL, "sigma")),
            targetData = c(1, 2, 3, 4)),
    "no column for"
  )

  ## a dataset of the wrong length
  expect_error(
    discFun(matrix(0, ncol = 1, dimnames = list(NULL, "mu")),
            targetData = c(1, 2)),
    "one value per data node"
  )
})
