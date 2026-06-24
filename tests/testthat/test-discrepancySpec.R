library(testthat)
library(nimble)
library(cppp)

## Check if we can make a simple discrepancy object
test_that("discrepancy() creates a basic discrepancy spec", {
  disc <- discrepancy("data_only", "mean")

  expect_s3_class(disc, "cppp_discrepancy")
  expect_equal(disc$kind, "data_only")
  expect_equal(disc$name, "mean")
  expect_null(disc$dataNodes)
  expect_null(disc$modelNodes)
})

## allow either one discrepancy or a list of discrepancies.
test_that("standardizeDiscrepancies() accepts one spec or a list", {
  disc1 <- discrepancy("data_only", "mean")
  disc2 <- discrepancy("data_plus_model", "deviance")

  out1 <- standardizeDiscrepancies(disc1)
  out2 <- standardizeDiscrepancies(list(disc1, disc2))

  expect_length(out1, 1)
  expect_length(out2, 2)
  expect_true(all(vapply(out2, inherits, logical(1), "cppp_discrepancy")))
})

## Here the model fills in the data nodes when the user leaves them out.
test_that("completeDiscrepancy() fills in default data nodes", {
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
  out <- completeDiscrepancy(model, disc)

  expect_equal(out$kind, "data_only")
  expect_equal(out$name, "mean")
  expect_equal(out$dataNodes, model$expandNodeNames("y", returnScalarComponents = TRUE))
})

## Some discrepancies need extra model-based quantities, so we check that.
test_that("completeDiscrepancy() requires model nodes when needed", {
  code <- nimbleCode({
    for (i in 1:n) {
      y[i] ~ dpois(lambda)
      y_exp[i] <- lambda
    }
    lambda ~ dgamma(1, 1)
  })

  model <- nimbleModel(
    code = code,
    constants = list(n = 3),
    data = list(y = c(1, 0, 2)),
    inits = list(lambda = 1)
  )

  expect_error(
    completeDiscrepancy(model, discrepancy("data_plus_model", "chisquared")),
    "requires `modelNodes`"
  )

  out <- completeDiscrepancy(
    model,
    discrepancy("data_plus_model", "chisquared", modelNodes = "y_exp")
  )

  expect_equal(out$modelNodes, model$expandNodeNames("y_exp", returnScalarComponents = TRUE))
})

## The NIMBLE discrepancy is the single source of truth. It runs uncompiled
## (R-interpreted) for development and compiled for speed; we check both against
## independent, hand-computed expected values.
test_that("makeDiscrepancyNimbleFun() evaluates simple discrepancies", {
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
  cModel <- compileNimble(model)

  ## Independent oracle: y = c(1,2,3,4), mu = 0, sd = 1.
  expectedMean <- 2.5
  expectedDev  <- -2 * sum(dnorm(c(1, 2, 3, 4), mean = 0, sd = 1, log = TRUE))

  meanNim <- makeDiscrepancyNimbleFun(model, discrepancy("data_only", "mean"))
  cMeanNim <- compileNimble(meanNim, project = cModel)

  devNim <- makeDiscrepancyNimbleFun(model, discrepancy("data_plus_model", "deviance"))
  cDevNim <- compileNimble(devNim, project = cModel)

  ## uncompiled (R-interpreted) execution
  expect_equal(meanNim$run(), expectedMean)
  expect_equal(devNim$run(), expectedDev)

  ## compiled execution agrees with uncompiled and with the oracle
  expect_equal(cMeanNim$run(), expectedMean)
  expect_equal(cDevNim$run(), expectedDev)
})
