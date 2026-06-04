library(testthat)
library(nimble)
library(cppp)

test_that("discrepancy() creates a basic discrepancy spec", {
  disc <- discrepancy("data_only", "mean")

  expect_s3_class(disc, "cppp_discrepancy")
  expect_equal(disc$kind, "data_only")
  expect_equal(disc$name, "mean")
  expect_null(disc$dataNodes)
  expect_null(disc$modelNodes)
  expect_null(disc$external)
})

test_that("standardizeDiscrepancies() accepts one spec or a list", {
  disc1 <- discrepancy("data_only", "mean")
  disc2 <- discrepancy("data_plus_model", "deviance")

  out1 <- standardizeDiscrepancies(disc1)
  out2 <- standardizeDiscrepancies(list(disc1, disc2))

  expect_length(out1, 1)
  expect_length(out2, 2)
  expect_true(all(vapply(out2, inherits, logical(1), "cppp_discrepancy")))
})

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
