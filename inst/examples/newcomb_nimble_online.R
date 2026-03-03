############################################################
## Newcomb data example: NIMBLE + online discrepancy + CPPP
############################################################
# Install a dev version alongside nimble using another folder
# library(devtools)

## If I want to install alomgside nimble I can use a specific path
## and load from there

# install_github("nimble-dev/nimble",
# 	subdir = "packages/nimble", #subdir
# 	ref = "derived_discrepancy",   #branch
# 	lib = "/opt/homebrew/Cellar/r/4.5.2_1/dev") #path to folder

## 1) Packages
library(nimble, lib.loc="/opt/homebrew/Cellar/r/4.5.2_1/dev")
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

## 4) online discrepancy
## Discrepancy: data-only, min(y)
## written using a nimble function

nfDef <- nimbleFunction(
  setup = TRUE,
  run = function(y = double(1)) {
    diff <- min(y)
    returnType(double())
    return(diff)
  }
)

Rnf <- nfDef()
conf <- configureMCMC(newcomb_model, monitors = paramNames)
conf$addDerivedQuantity(derived_discrepancy,
                        control = list(simNodes = 'y', discrepancyFunction = Rnf))
conf$printDerivedQuantities()

Rmcmc <- buildMCMC(conf)

## short run to look into the output
# test <- runMCMC(Rmcmc, niter = 20)
# test$derived

## 5) define a function that generates new data

control <- list(
  model      = newcomb_model,
  dataNames  = dataNames,
  paramNames = paramNames,
  verbose    = TRUE
)


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
# Compile and run MCMC

Cnewcomb_model <- compileNimble(newcomb_model)
Cmcmc <- compileNimble(Rmcmc, project = newcomb_model)

out <- runMCMC(
  Cmcmc,
  niter   = 5000,
  nburnin = 1000
)
is(out)

MCMCSamples <- as.matrix(out$samples)
discSamples <- as.matrix(out$derived$discrepancy)
head(MCMCSamples)
head(discSamples)

