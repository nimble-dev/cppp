library(testthat)
library(nimble)
library(cppp)

## Check if we can make a simple discrepancy object
test_that("discrepancy() creates a basic discrepancy spec", {
  disc <- discrepancy("mean")

  expect_s3_class(disc, "cppp_discrepancy")
  expect_equal(disc$name, "mean")
  expect_null(disc$dataNodes)
  expect_null(disc$modelNodes)
  expect_null(disc$fun)
})

## Without `fun`, the name has to be one the package knows, and the error should
## say what the options are.
test_that("discrepancy() rejects an unknown name with no `fun`", {
  expect_error(discrepancy("notADiscrepancy"), "Unknown discrepancy")
  expect_error(discrepancy("notADiscrepancy"), "mean")
})

## `name` labels the results, so a user-supplied discrepancy must not borrow a
## built-in's name and compute something else under it.
test_that("discrepancy() accepts a user-supplied fun but not a built-in name", {
  myDisc <- nimbleFunction(
    setup = function(model, dataNodes, modelNodes) {},
    run = function() {
      returnType(double(0))
      return(min(values(model, dataNodes)))
    }
  )

  disc <- discrepancy("mymin", fun = myDisc)
  expect_s3_class(disc, "cppp_discrepancy")
  expect_equal(disc$name, "mymin")
  expect_false(is.null(disc$fun))

  expect_error(discrepancy("mean", fun = myDisc), "built-in discrepancy")
  expect_error(discrepancy("mymin", fun = "not a function"), "must be a nimbleFunction")
})

## allow either one discrepancy or a list of discrepancies.
test_that("standardizeDiscrepancies() accepts one spec or a list", {
  disc1 <- discrepancy("mean")
  disc2 <- discrepancy("deviance")

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

  disc <- discrepancy("mean")
  out <- completeDiscrepancy(model, disc)

  expect_equal(out$name, "mean")
  expect_equal(out$dataNodes, model$expandNodeNames("y", returnScalarComponents = TRUE))
  ## a discrepancy that does not use modelNodes still gets them, as character(0)
  expect_equal(out$modelNodes, character(0))
})

## Some discrepancies need extra model-based quantities. The requirement lives
## with the discrepancy, so it is reported when the discrepancy is built.
test_that("a discrepancy that needs model nodes complains when they are missing", {
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
    makeDiscrepancyNimbleFun(model, discrepancy("chisquared")),
    "requires `modelNodes`"
  )

  out <- completeDiscrepancy(
    model,
    discrepancy("chisquared", modelNodes = "y_exp")
  )

  expect_equal(out$modelNodes, model$expandNodeNames("y_exp", returnScalarComponents = TRUE))
})

## The NIMBLE discrepancy runs uncompiled
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

  meanNim <- makeDiscrepancyNimbleFun(model, discrepancy("mean"))
  cMeanNim <- compileNimble(meanNim, project = cModel)

  devNim <- makeDiscrepancyNimbleFun(model, discrepancy("deviance"))
  cDevNim <- compileNimble(devNim, project = cModel)

  ## uncompiled (R-interpreted) execution
  expect_equal(meanNim$run(), expectedMean)
  expect_equal(devNim$run(), expectedDev)

  ## compiled execution agrees with uncompiled and with the oracle
  expect_equal(cMeanNim$run(), expectedMean)
  expect_equal(cDevNim$run(), expectedDev)
})

## A user-supplied discrepancy must be built and run exactly like a built-in,
## uncompiled and compiled.
test_that("makeDiscrepancyNimbleFun() builds a user-supplied discrepancy", {
  code <- nimbleCode({
    for (i in 1:n) {
      y[i] ~ dnorm(mu, sd = 1)
    }
    mu ~ dnorm(0, sd = 10)
  })

  model <- nimbleModel(
    code = code,
    constants = list(n = 4),
    data = list(y = c(1, 2, 3, 10)),
    inits = list(mu = 2)
  )
  cModel <- compileNimble(model)

  ## uses modelNodes, to check they are passed through
  myGap <- nimbleFunction(
    setup = function(model, dataNodes, modelNodes) {},
    run = function() {
      y  <- values(model, dataNodes)
      mu <- values(model, modelNodes)[1]
      returnType(double(0))
      return(max(y) - mu)
    }
  )

  disc <- discrepancy("mygap", modelNodes = "mu", fun = myGap)
  nim <- makeDiscrepancyNimbleFun(model, disc)
  cNim <- compileNimble(nim, project = cModel)

  expect_equal(nim$run(), 8)   # max(1,2,3,10) - 2
  expect_equal(cNim$run(), 8)
})
