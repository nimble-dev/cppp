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
  expect_null(disc$external)
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

## Now we make a plain R version and see it run.
test_that("makeDiscrepancyRFun() evaluates simple discrepancies", {
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

  meanFun <- makeDiscrepancyRFun(model, discrepancy("data_only", "mean"))
  devFun <- makeDiscrepancyRFun(model, discrepancy("data_plus_model", "deviance"))

  expect_equal(meanFun(), 2.5)
  expect_true(is.numeric(devFun()))
  expect_length(devFun(), 1)

  expect_equal(meanFun(c(2, 2, 2, 2)), 2)
})

## Finally, we make a NIMBLE version from the same discrepancy and check
## that it gives the same answer as the plain R version.
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

  meanR <- makeDiscrepancyRFun(model, discrepancy("data_only", "mean"))
  meanNim <- makeDiscrepancyNimbleFun(model, discrepancy("data_only", "mean"))
  cMeanNim <- compileNimble(meanNim, project = cModel)

  devR <- makeDiscrepancyRFun(model, discrepancy("data_plus_model", "deviance"))
  devNim <- makeDiscrepancyNimbleFun(model, discrepancy("data_plus_model", "deviance"))
  cDevNim <- compileNimble(devNim, project = cModel)

  expect_equal(cMeanNim$run(), meanR())
  expect_equal(cDevNim$run(), devR())
})
