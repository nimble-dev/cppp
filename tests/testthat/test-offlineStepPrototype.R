library(testthat)
library(nimble)
library(cppp)

## This test shows one full step:
## set a posterior-draw state, compute the observed discrepancy,
## simulate a replicate, then compute the replicated discrepancy.
test_that("runOneDiscrepancyStepR() runs one discrepancy step", {
  code <- nimbleCode({
    for (i in 1:n) {
      y[i] ~ dnorm(mu, sd = 1)
    }
    mu ~ dnorm(0, sd = 10)
  })

  model <- nimbleModel(
    code = code,
    constants = list(n = 4),
    data = list(y = c(1, 2, 3, 4)),
    inits = list(mu = 0)
  )

  disc <- discrepancy("data_only", "mean")
  sim <- simulation("conditional")

  out <- runOneDiscrepancyStepR(
    model = model,
    disc = disc,
    sim = sim,
    thetaRow = c(mu = 1.5),
    paramNodes = "mu"
  )

  expect_true(is.list(out))
  expect_true(all(c("obs", "rep") %in% names(out)))
  expect_true(is.numeric(out$obs))
  expect_true(is.numeric(out$rep))
  expect_length(out$obs, 1)
  expect_length(out$rep, 1)
})
