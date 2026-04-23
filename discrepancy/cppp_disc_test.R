library(nimble)

# read in disc functions
source("discrepancy/cppp_discrepancy.R")
source("discrepancy/cppp_discrepancy_NIMBLE.R")

####################################
# Set up model, data, MCMC samples #
####################################

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


#####################
# Get discrepancies #
#####################

# choose discType
discType <- c("mean", "deviance", "chisq")

# disc control (updated from main code)
discControl = list(
  model         = model,
  dataNames     = dataNames,
  dataNodes     = dataNodes,
  paramNames    = paramNames,
  paramNodes    = paramNodes,
  expectedNames = expectedNames, # new
  expectedNodes = expectedNodes, # new
  discType      = discType # new
)

####
# Non-NIMBLE functions 
####

# calculate discrepancies
obsDisc <- calcDiscFunction(MCMCSamples  = MCMCSamples,
                            targetData   = observedData,
                            control      = discControl)

# scalar PPP for the observed data
obsPPP <- sapply(obsDisc, function(x) {
  mean(as.numeric(x[, "obs"]) >= as.numeric(x[, "sim"]))
})


####
# NIMBLE functions 
####

# create instance of discrepancy function
calc_disc <- calcDiscFunction_nimble(MCMCSamples, as.numeric(observedData),
                                     discControl)

# compile
Ccalc_disc <- compileNimble(calc_disc, project = model_uncompiled)

# run - this keeps crashing R :/
obsDisc_nimble <- Ccalc_disc$run()
