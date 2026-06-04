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

## In marginal mode, the user must currently say which nodes to resimulate.
test_that("completeSimulation() requires simulate nodes in marginal mode", {
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
    "simulateNodes"
  )

  sim <- simulation(mode = "marginal", simulateNodes = "y")
  out <- completeSimulation(model, sim)

  expect_equal(out$mode, "marginal")
  expect_equal(out$dataNodes, model$expandNodeNames("y", returnScalarComponents = TRUE))
  expect_equal(out$simulateNodes, model$expandNodeNames("y", returnScalarComponents = TRUE))
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
