library(testthat)
library(nimble)
library(cppp)

## Marginal mode, end to end: the latents b are redrawn for each replicate.

makeMarginalTestModel <- function(n = 200) {
  code <- nimbleCode({
    for (i in 1:n) {
      b[i] ~ dnorm(0, sd = tau)
      y[i] ~ dnorm(mu + b[i], sd = 1)
    }
    mu  ~ dnorm(0, sd = 10)
    tau ~ dunif(0, 20)
  })
  nimbleModel(code, constants = list(n = n), data = list(y = rep(0, n)),
              inits = list(mu = 0, tau = 1, b = rep(0, n)))
}

## With tau = 10, a replicate that redraws b has spread about sqrt(1 + 10^2).
## One that keeps b fixed has spread about 1.
test_that("makeSimulateNewDataFun() redraws the latents in marginal mode", {
  model <- makeMarginalTestModel()

  marg <- makeSimulateNewDataFun(
    model$newModel(),
    completeNodes(model, paramNodes = c("mu", "tau"),
                  simulation = simulation("marginal", latentNodes = "b"))
  )
  cond <- makeSimulateNewDataFun(
    model$newModel(),
    completeNodes(model, paramNodes = c("mu", "tau"))
  )

  set.seed(1)
  expect_gt(sd(marg(c(mu = 0, tau = 10))), 5)
  expect_lt(sd(cond(c(mu = 0, tau = 10))), 2)
})

test_that("runCalibrationNIMBLE() runs in marginal mode", {
  model <- makeMarginalTestModel(n = 20)

  draws <- cbind(mu  = c(0, 0.5, 1),
                 tau = c(1, 1.5, 2))

  set.seed(1)
  res <- runCalibrationNIMBLE(
    model          = model,
    dataNames      = "y",
    paramNames     = c("mu", "tau"),
    MCMCSamples    = draws,
    discrepancies  = list(discrepancy("mean"), discrepancy("deviance")),
    simulation     = simulation("marginal", latentNodes = "b"),
    nReps          = 2,
    MCMCcontrolRep = list(niter = 10, nburnin = 0, thin = 1)
  )

  expect_named(res$CPPP, c("mean", "deviance"))
  expect_named(res$obsPPP, c("mean", "deviance"))
  expect_equal(dim(res$repPPP), c(2L, 2L))
  expect_true(all(is.finite(res$repPPP)))
})
