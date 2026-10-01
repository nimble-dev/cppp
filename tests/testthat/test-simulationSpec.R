library(testthat)
library(nimble)
library(cppp)

## Check if we can make a simple simulation object
test_that("simulation() creates a basic simulation spec", {
  sim <- simulation()

  expect_s3_class(sim, "cppp_simulation")
  expect_equal(sim$mode, "conditional")
  expect_null(sim$latentNodes)

  expect_error(simulation("conditional", latentNodes = "b"), "only used")
})

makeLatentModel <- function() {
  code <- nimbleCode({
    for (i in 1:n) {
      b[i] ~ dnorm(0, sd = tau)
      y[i] ~ dnorm(mu + b[i], sd = 1)
    }
    mu ~ dnorm(0, sd = 10)
    tau ~ dunif(0, 10)
  })
  nimbleModel(code, constants = list(n = 3), data = list(y = c(1, 2, 3)),
              inits = list(mu = 0, tau = 1, b = rep(0, 3)))
}

## In marginal mode only the user knows the latent states; there is no default.
test_that("completeNodes() requires latent nodes in marginal mode", {
  model <- makeLatentModel()

  expect_error(completeNodes(model, simulation = simulation("marginal")),
               "latentNodes")
  expect_error(completeNodes(model, simulation = simulation("marginal", latentNodes = "y")),
               "must not include data")
})

## The difference between the two modes, on a model with latent effects.
test_that("marginal redraws the latents and everything below, conditional only the data", {
  model <- makeLatentModel()
  bNodes <- model$expandNodeNames("b", returnScalarComponents = TRUE)
  yNodes <- model$expandNodeNames("y", returnScalarComponents = TRUE)

  marg <- completeNodes(model, simulation = simulation("marginal", latentNodes = "b"))
  cond <- completeNodes(model, simulation = simulation("conditional"))

  ## the user names b; the helper nodes NIMBLE puts between b and y, and y
  ## itself, are added
  expect_true(all(bNodes %in% marg$simulate))
  expect_true(all(yNodes %in% marg$simulate))
  expect_true(any(grepl("^lifted_", marg$simulate)))

  expect_false(any(bNodes %in% cond$simulate))
  expect_equal(cond$simulate, cond$data)
})

## Naming the top latent layer redraws the lower ones too, past the first
## stochastic node below it.
test_that("completeNodes() redraws every layer below the named latents", {
  code <- nimbleCode({
    for (i in 1:n) {
      b[i] ~ dnorm(0, sd = 1)
      c[i] ~ dnorm(b[i], sd = 1)
      y[i] ~ dnorm(c[i], sd = 1)
    }
  })
  model <- nimbleModel(code, constants = list(n = 2), data = list(y = c(1, 2)),
                       inits = list(b = c(0, 0), c = c(0, 0)))

  out <- completeNodes(model, paramNodes = "b",
                       simulation = simulation("marginal", latentNodes = "b"))

  expect_setequal(out$simulate, c("b[1]", "b[2]", "c[1]", "c[2]", "y[1]", "y[2]"))
})

test_that("completeNodes() fills in the data and parameter defaults", {
  model <- makeLatentModel()
  out <- completeNodes(model)

  expect_s3_class(out, "cppp_nodes")
  expect_equal(out$data, model$expandNodeNames("y", returnScalarComponents = TRUE))
  expect_setequal(out$params, model$expandNodeNames(c("b", "mu", "tau"),
                                                    returnScalarComponents = TRUE))
  expect_true(all(c(out$params, out$paramDeps, out$simulate) %in% out$saved))

  expect_error(completeNodes(model, dataNodes = "b", paramNodes = character(0)),
               "did not match")
})
