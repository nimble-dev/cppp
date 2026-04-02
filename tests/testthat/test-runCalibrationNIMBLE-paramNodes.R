library(testthat)
library(nimble)
library(cppp)

### model setup - same for all the tests
set.seed(1)

n <- 20
x <- rnorm(n)
beta_true <- c(1, -0.5)
y <- rnorm(n, mean = beta_true[1] + beta_true[2] * x, sd = 1)

regCode <- nimbleCode({
  for (i in 1:n) {
    mu[i] <- beta[1] + beta[2] * x[i]
    y[i] ~ dnorm(mu[i], sd = 1)
  }
  for (j in 1:2) {
    beta[j] ~ dnorm(0, sd = 10)
  }
})

baseModel <- nimbleModel(
  code = regCode,
  constants = list(n = n, x = x),
  data = list(y = y),
  inits = list(beta = c(0, 0))
)

dataNames  <- "y"
paramNames <- "beta"
paramNodes <- baseModel$expandNodeNames(paramNames)

simDataFun <- function(thetaRow, control) {
  model      <- control$model
  dataNodes  <- control$dataNodes
  paramNodes <- control$paramNodes

  values(model, paramNodes) <- thetaRow[paramNodes]
  model$simulate(nodes = dataNodes, includeData = TRUE)
  values(model, dataNodes)
}

mean_disc <- function(data, thetaRow, ...) {
  mean(data)
}

discConfig <- list(
  simulateNewDataFun = simDataFun,
  discrepancy = mean_disc
)

discFun <- makeOfflineDiscFun(discConfig)

test_that("runCalibrationNIMBLE works with expanded paramNames in main MCMC", {
  model <- baseModel$newModel()

  control <- list(
    model   = model,
    verbose = FALSE
  )

  res <- runCalibrationNIMBLE(
    model = model,
    dataNames = dataNames,
    paramNames = paramNames,
    discFun = discFun,
    simulateNewDataFun = simDataFun,
    nReps = 2,
    MCMCcontrolMain = list(niter = 100, nburnin = 20, thin = 1),
    MCMCcontrolRep  = list(niter = 30, nburnin = 0, thin = 1),
    control = control
  )

  expect_true(is.list(res))
  expect_true("CPPP" %in% names(res))
  expect_true("repPPP" %in% names(res))
})

test_that("runCalibrationNIMBLE works with supplied MCMCSamples and expanded paramNames", {
  model_for_mcmc <- baseModel$newModel()

  samples <- nimbleMCMC(
    model_for_mcmc,
    niter = 120,
    nburnin = 20,
    monitors = paramNodes
  )
  MCMCSamples <- as.matrix(samples)

  model <- baseModel$newModel()

  control <- list(
    model   = model,
    verbose = FALSE
  )

  res <- runCalibrationNIMBLE(
    model = model,
    dataNames = dataNames,
    paramNames = paramNames,
    MCMCSamples = MCMCSamples,
    discFun = discFun,
    simulateNewDataFun = simDataFun,
    nReps = 2,
    MCMCcontrolRep  = list(niter = 30, nburnin = 0, thin = 1),
    control = control
  )

  expect_true(is.list(res))
  expect_true(all(paramNodes %in% colnames(MCMCSamples)))
  expect_true("CPPP" %in% names(res))
  expect_true("repPPP" %in% names(res))
})

test_that("runCalibrationNIMBLE errors when expanded param nodes are missing from MCMCSamples", {
  model <- baseModel$newModel()

  badSamples <- cbind(`beta[1]` = rnorm(20))
  badSamples <- as.matrix(badSamples)

  control <- list(
    model   = model,
    verbose = FALSE
  )

  expect_error(
    runCalibrationNIMBLE(
      model = model,
      dataNames = dataNames,
      paramNames = paramNames,
      MCMCSamples = badSamples,
      discFun = discFun,
      simulateNewDataFun = simDataFun,
      nReps = 2,
      control = control
    ),
    "paramNames missing from MCMCSamples"
  )
})
