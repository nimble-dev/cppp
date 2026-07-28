library(testthat)
library(nimble)
library(cppp)

## Check if we can make a simple simulation object
test_that("simulation() creates a basic simulation spec", {
  sim <- simulation()

  expect_s3_class(sim, "cppp_simulation")
  expect_equal(sim$mode, "conditional")
  expect_null(sim$dataNodes)
  expect_null(sim$simulateNodes)
})

## In marginal mode we need to know the parameters, unless the user has already
## said which nodes to resimulate.
test_that("completeSimulation() needs paramNodes or simulateNodes in marginal mode", {
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

  expect_error(
    completeSimulation(model, simulation(mode = "marginal")),
    "paramNodes"
  )

  sim <- simulation(mode = "marginal", simulateNodes = "y")
  out <- completeSimulation(model, sim)

  expect_equal(out$mode, "marginal")
  expect_equal(out$dataNodes, model$expandNodeNames("y", returnScalarComponents = TRUE))
  expect_equal(out$simulateNodes, model$expandNodeNames("y", returnScalarComponents = TRUE))
})

## The point of marginal mode: latent nodes are redrawn too, not just the data.
test_that("completeSimulation() works out marginal simulate nodes from paramNodes", {
  code <- nimbleCode({
    for (i in 1:n) {
      b[i] ~ dnorm(0, sd = tau)
      y[i] ~ dnorm(mu + b[i], sd = 1)
    }
    mu ~ dnorm(0, sd = 10)
    tau ~ dunif(0, 10)
  })

  model <- nimbleModel(
    code = code,
    constants = list(n = 3),
    data = list(y = c(1, 2, 3)),
    inits = list(mu = 0, tau = 1, b = rep(0, 3))
  )

  paramNodes <- c("mu", "tau")

  marg <- completeSimulation(model, simulation(mode = "marginal"), paramNodes)
  cond <- completeSimulation(model, simulation(mode = "conditional"), paramNodes)

  expected <- model$expandNodeNames(c("b", "y"), returnScalarComponents = TRUE)
  expect_setequal(marg$simulateNodes, expected)

  ## conditional keeps the latent effects fixed; marginal does not
  expect_setequal(cond$simulateNodes, cond$dataNodes)
  expect_true(all(model$expandNodeNames("b", returnScalarComponents = TRUE) %in%
                    marg$simulateNodes))
  expect_false(any(model$expandNodeNames("b", returnScalarComponents = TRUE) %in%
                     cond$simulateNodes))

  ## parents come before children, so simulating in this order is safe
  expect_lt(
    max(match(model$expandNodeNames("b", returnScalarComponents = TRUE), marg$simulateNodes)),
    min(match(model$expandNodeNames("y", returnScalarComponents = TRUE), marg$simulateNodes))
  )
})

## In conditional mode, the default is to simulate only the data nodes.
test_that("completeSimulation() uses data nodes as default simulate nodes in conditional mode", {
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

  sim <- simulation(mode = "conditional")
  out <- completeSimulation(model, sim)

  expect_equal(out$simulateNodes, out$dataNodes)
})
