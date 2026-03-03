# tests/testthat/test-runCalibration.R

test_that("runCalibration returns numeric PPP vector of correct length", {
  skip_if_not("runCalibration" %in% ls("package:cpppPrototype"),
              "runCalibration() not exported yet")

  set.seed(1)

  observedData <- rnorm(20)
  # Fake "long-run" MCMC samples
  MCMCSamples <- matrix(0, nrow = 10, ncol = 1)

  # Toy simulateNewDataFun: ignore MCMCSamples, just simulate data
  simulateNewDataFun <- function(MCMCSamples, ...) {
    rnorm(20)
  }

  # Toy MCMCFun: given data, return a "chain" whose mean depends on the data
  MCMCFun <- function(newData, control, ...) {
    theta <- rnorm(100, mean(newData), sd = 1 / sqrt(length(newData)))
    cbind(theta = theta)
  }

  # Toy discFun: discrepancy is |mean(y)|; PPP via comparing |theta| to |mean(y)|
  discFun <- function(MCMCSamples, ...) {
    discObs <- abs(mean(MCMCSamples[, "theta"]))
    discRep <- abs(MCMCSamples[, "theta"])  # pretend these are replicated discrepancies
    list(
      rep = discRep,
      obs = discObs
    )
  }

  nReps <- 5

  repPPP <- runCalibration(
    MCMCSamples = MCMCSamples,
    observedData = observedData,
    MCMCFun     = MCMCFun,
    simulateNewDataFun = simulateNewDataFun,
    discFun     = discFun,
    nReps     = nReps
  )

  expect_type(repPPP, "double")
  expect_length(repPPP, nReps)
  expect_true(all(is.finite(repPPP)))
  expect_true(all(repPPP >= 0 & repPPP <= 1))
})

test_that("runCalibration checks nReps and discFun output", {
  skip_if_not("runCalibration" %in% ls("package:cpppPrototype"),
              "runCalibration() not exported yet")

  observedData <- rnorm(20)
  MCMCSamples <- matrix(0, nrow = 5, ncol = 1)

  simulateNewDataFun <- function(MCMCSamples, ...) rnorm(10)
  MCMCFun     <- function(newData, control, ...) matrix(newData, ncol = 1)

  bad_discFun <- function(MCMCSamples, ...) {
    0.5  # numeric scalar, no $rep / $obs
  }

  # nReps must be positive integer
  expect_error(
    runCalibration(MCMCSamples,
                   observedData,
                   MCMCFun,
                   simulateNewDataFun,
                   bad_discFun,
                   nReps = 0),
    "positive integer"
  )

  # discFun must return list with $rep and $obs
  expect_error(
    runCalibration(MCMCSamples,
                   observedData, MCMCFun,
                   simulateNewDataFun,
                   bad_discFun,
                   nReps = 1),
    "list with components `rep` and `obs`"
  )
})
