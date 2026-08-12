
########################################
# create sample model and MCMC samples #
########################################

# occupancy model
occu_model <- nimbleCode({

  psi ~ dbeta(1, 1)
  p ~ dbeta(1, 1)

  for (i in 1:nSites) {

    z[i] ~ dbern(psi)

    # observed data
    y[i] ~ dbinom(z[i] * p, nVisits)

    # expected data
    y_exp[i] <- z[i] * p * nVisits

  }

})

# function for simulating data
simulate_occu <- function(params, nSites, nVisits) {

  y <- rep(NA, nSites)

  z <- rbinom(nSites, 1, params$psi)

  for (i in 1:nSites) {

    y[i] <- rbinom(1, nVisits, z[i] * params$p)

  }

  return(y)
}

# model components
nVisits <- 6
nSites <- 20
params <- list(p = 0.3, psi = 0.4)
constants <- list(nSites = nSites, nVisits = nVisits)
observedData <- simulate_occu(params = list(p = 0.3, psi = 0.4),
                              nSites = 20, nVisits = 6)

# uncompiled model
model_uncompiled <- nimbleModel(occu_model, constants = constants,
                                data = list(y = observedData))

# compiled model
model <- compileNimble(model_uncompiled)

# get data, param, and expected nodes
dataNames <- "y"
dataNodes <- model$expandNodeNames(dataNames)

paramNames <- c("p", "psi", "z")
paramNodes <- model$expandNodeNames(paramNames)

expectedNames <- "y_exp"
expectedNodes <- model$expandNodeNames(expectedNames)

# get MCMC samples
# MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1)
MCMCcontrolMain = list(niter = 500, nburnin = 100, thin = 1)

mcmcConfFun <- function(model) {
  configureMCMC(model, monitors = paramNodes, print = FALSE)
}
mcmcConf       <- mcmcConfFun(model)
mcmcUncompiled <- buildMCMC(mcmcConf)
cmcmc          <- compileNimble(mcmcUncompiled, project = model,
                                resetFunctions = TRUE)

obsMCMC <- runMCMC(
  cmcmc,
  niter   = MCMCcontrolMain$niter,
  nburnin = MCMCcontrolMain$nburnin,
  thin    = MCMCcontrolMain$thin
)
MCMCSamples <- as.matrix(obsMCMC)

# test discrepancy calculation output
test_that("GenericDiscrFunction works for basic cases", {

  # get original data and param values
  original_data <- values(model, dataNodes)
  original_param <- values(model, paramNodes)

  # choose discType
  discType <- c("mean", "variance", "deviance", "chisquared", "freemantukey")

  # disc control
  discControl = list(
    model         = model_uncompiled,
    dataNames     = dataNames,
    dataNodes     = dataNodes,
    paramNames    = paramNames,
    paramNodes    = paramNodes,
    expectedNames = expectedNames,
    expectedNodes = expectedNodes,
    discType      = discType
  )

  # create instance of discrepancy function
  calc_disc <- calcDiscrepancies(discControl)

  # compile
  Ccalc_disc <- compileNimble(calc_disc, project = model_uncompiled)

  # run
  out <- Ccalc_disc$run(MCMCSamples)

  # test dimensions of output
  expect_equal(dim(out), c(5, 400, 2))

  # test that no values are NA
  expect_true(all(!is.na(out)))

  # test that original data values have been returned to model
  expect_equal(model$y, original_data)
  expect_equal(model_uncompiled$y, original_data)

  # test that original parameter values have been returned to model
  expect_equal(values(model, paramNames), original_param)

})

# test error trapping
test_that("GenericDiscrFunction error trap works", {

  # choose discType
  discType <- c("mean", "deviance", "chisq")

  # disc control
  discControl = list(
    model         = model_uncompiled,
    dataNames     = dataNames,
    dataNodes     = dataNodes,
    paramNames    = paramNames,
    paramNodes    = paramNodes,
    expectedNames = expectedNames,
    expectedNodes = expectedNodes,
    discType      = discType
  )

  # create instance of discrepancy function
  expect_error(
    calcDiscrepancies(discControl),
    paste0("discType 'chisq' is invalid. Must be: mean, variance, deviance, ",
           "chisquared, or freemantukey")
  )

})
