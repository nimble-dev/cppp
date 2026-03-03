############################################################
## Newcomb data example: NIMBLE + offline discrepancy + CPPP
############################################################

## 1) Packages
library(nimble)
library(cppp)

## 2) Data: Newcomb light-speed measurements
lightPath <- system.file("examples", "light.txt", package = "cppp")
newcombData <- list(y = read.table(lightPath)$V1)

constants <- list(n = length(newcombData$y))


## 3) NIMBLE model
newcombCode <- nimbleCode({
  for (i in 1:n) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
  mu ~ dflat()
  log(sigma) ~ dflat()
})


inits <- list(mu = 0, log_sigma = 2)

newcomb_model <- nimbleModel(
  code      = newcombCode,
  data      = newcombData,
  inits     = inits,
  constants = constants
)

dataNames  <- "y"
paramNames <- c("mu", "sigma")

## 4) Offline discrepancy
## Discrepancy: data-only, min(y)

min_disc <- function(data, thetaRow, control) {
  min(data)
}


## asymmetry discrepancy using R
asymm_disc <- function(data, thetaRow, control) {
  mu <- thetaRow["mu"]
  dataSorted <- sort(data)
  abs(dataSorted[61] - mu) - abs(dataSorted[6] - mu)
}

control <- list(
  model      = newcomb_model,
  dataNames  = dataNames,
  paramNames = paramNames,
  verbose    = TRUE
)

## function that generates new data
newcombNewData <- function(thetaRow, control) {
  model      <- control$model
  dataNames  <- control$dataNames
  paramNames <- control$paramNames

  for (nm in paramNames) {
    model[[nm]] <- thetaRow[nm]
  }

  model$simulate(nodes = dataNames, includeData = TRUE)
  model[[dataNames]]
}

##########################################
# 0) Run MCMC
samples <- nimbleMCMC(
  newcomb_model,
  niter   = 5000,
  nburnin = 1000,
  monitors = paramNames
)
str(samples)
MCMCSamples <- as.matrix(samples)
head(MCMCSamples)

##  Check simulation
# control <- list(
#   model      = newcomb_model,
#   dataNames  = dataNames,
#   paramNames = paramNames
# )
#
# theta1 <- MCMCSamples[1, , drop = TRUE]
# ySim1 <- newcombNewData(thetaRow = theta1, control = control)
#
# str(ySim1)
# length(ySim1)
########################################################
## Build disc_fun via the package helper

discConfig <- list(
  simulateNewDataFun = newcombNewData,
  discrepancy        = asymm_disc
)

discFun <- makeOfflineDiscFun(discConfig)
set.seed(1)

resNewcomb <- runCalibrationNIMBLE(
  model = newcomb_model$newModel(),
  dataNames = dataNames,
  paramNames = paramNames,
  discFun = discFun,
  simulateNewDataFun = newcombNewData,
  nReps = 10,
  MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1),
  MCMCcontrolRep  = list(niter = 100,  nburnin = 0,    thin = 1),
  control = control
)

print(resNewcomb$CPPP)
print(resNewcomb$repPPP)


#####################
### check providing the MCMC chain first

resNewcomb <- runCalibrationNIMBLE(
  model = newcomb_model$newModel(),
  dataNames = dataNames,
  paramNames = paramNames,
  MCMCSamples = MCMCSamples,
  discFun = discFun,
  simulateNewDataFun = newcombNewData,
  nReps = 10,
  MCMCcontrolRep  = list(niter = 100,  nburnin = 0,    thin = 1),
  control = control
)

print(resNewcomb$CPPP)
print(resNewcomb$repPPP)

# obsDisc <- res_newcomb$discrepancies$obs
# plot(obsDisc$obs, obsDisc$sim)
# abline()
############################

#
# ##  Check discFun on the original (real) data
# disc <- discFun(
#   MCMCSamples = MCMCSamples,
#   targetData  = newcombData$y,
#   control     = control
# )
# plot(disc$obs, disc$sim,
#      xlab = "D(y, θ)",
#      ylab = "D(y*, θ)")
# abline(0, 1)
#
# mean(disc$sim >= disc$obs)
