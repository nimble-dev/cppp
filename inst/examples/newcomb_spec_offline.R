############################################################
## Newcomb data example: discrepancy specifications + CPPP
############################################################
##
## The discrepancies are written as specifications and turned into a
## calculator by the package, instead of being hand-written R functions.
##
## Two discrepancies at once:
##   mean   -- a built-in, data only
##   asymm  -- your own, and it needs mu

library(nimble)
library(cppp)

## 1) Data: Newcomb light-speed measurements
lightPath   <- system.file("examples", "light.txt", package = "cppp")
newcombData <- list(y = read.table(lightPath)$V1)
constants   <- list(n = length(newcombData$y))

## 2) Model
newcombCode <- nimbleCode({
  for (i in 1:n) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
  mu ~ dflat()
  log(sigma) ~ dflat()
})

newcombModel <- nimbleModel(
  code      = newcombCode,
  data      = newcombData,
  inits     = list(mu = 0, log_sigma = 2),
  constants = constants
)

dataNames  <- "y"
paramNames <- c("mu", "sigma")

## 3) Your own discrepancy: asymmetry of the two tails around mu.
##    Wrap R's sort() so the same code also works once compiled.
sortR <- nimbleRcall(prototype  = function(x = double(1)) {},
                     Rfun       = "sort",
                     returnType = double(1))

asymmDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {},
  run = function() {
    y  <- values(model, dataNodes)
    mu <- values(model, modelNodes)[1]
    ys <- sortR(y)
    returnType(double(0))
    ## 6th and 61st order statistics of the 66 measurements
    return(abs(ys[61] - mu) - abs(ys[6] - mu))
  }
)

discs <- list(
  discrepancy("mean"),
  discrepancy("asymm", modelNodes = "mu", fun = asymmDisc)
)

sim <- simulation("conditional")

## 4) Posterior draws for the real data
set.seed(1)
samples     <- nimbleMCMC(newcombModel, niter = 5000, nburnin = 1000,
                          monitors = paramNames)
MCMCSamples <- as.matrix(samples)

## Thinned only to keep this example quick;
thinned <- MCMCSamples[seq(1, nrow(MCMCSamples), by = 20), , drop = FALSE]

############################################################
## Stage 1: the calculator on its own, before any calibration
############################################################

## Its own copy of the model: it writes into the model as it works. This
## compiles the model and the two discrepancies, which takes a moment once.
calc <- makeDiscrepancyCalculator(
  model         = newcombModel$newModel(),
  discrepancies = discs,
  simulation    = sim,
  paramNodes    = paramNames
)

d <- calc(thinned, targetData = newcombData$y)

str(d)
## observed posterior predictive p-value, one per discrepancy
colMeans(d$sim >= d$obs)

## A normal model fits the middle of the data but not its left tail, so we
## expect mean near 0.5 and asymm far from it.

# par(mfrow = c(1, 2))
# for (nm in colnames(d$obs)) {
#   plot(d$obs[, nm], d$sim[, nm], xlab = "D(y, theta)", ylab = "D(y*, theta)",
#        main = nm)
#   abline(0, 1)
# }

############################################################
## Stage 2: cross-check the observed side against plain R
############################################################
## D(y, theta) has no randomness in it, so an independent implementation must
## give the same numbers.

asymmR <- function(data, thetaRow) {
  mu <- thetaRow["mu"]
  ys <- sort(data)
  abs(ys[61] - mu) - abs(ys[6] - mu)
}

obsAsymmR <- apply(thinned, 1, function(th) asymmR(newcombData$y, th))
max(abs(obsAsymmR - d$obs[, "asymm"]))     # expect ~0
max(abs(mean(newcombData$y) - d$obs[, "mean"]))  # expect 0

############################################################
## Stage 3: a small calibration
############################################################
## Hand the same specifications straight to runCalibrationNIMBLE. It builds
## both the calculator and the replicate-simulating function itself, so there
## is nothing model-specific left to write here.
##
## Note this compiles the discrepancies a second time: `calc` above is not
## reused, runCalibrationNIMBLE makes its own on its own copy of the model.

set.seed(1)
resNewcomb <- runCalibrationNIMBLE(
  model          = newcombModel$newModel(),
  dataNames      = dataNames,
  paramNames     = paramNames,
  MCMCSamples    = thinned,
  discrepancies  = discs,
  simulation     = sim,
  nReps          = 10,
  MCMCcontrolRep = list(niter = 10, nburnin = 0, thin = 1)
)

resNewcomb$obsPPP    # one per discrepancy
resNewcomb$CPPP      # one per discrepancy
resNewcomb$repPPP    # nReps x 2

#######
# nReps <- nrow(resNewcomb$repPPP)
# se    <- sqrt(resNewcomb$CPPP * (1 - resNewcomb$CPPP) / nReps)
#
# data.frame(
#   discrepancy = names(resNewcomb$CPPP),
#   obsPPP = round(resNewcomb$obsPPP, 4),
#   CPPP   = round(resNewcomb$CPPP, 4),
#   se     = round(se, 4),
#   lo95   = round(pmax(0, resNewcomb$CPPP - 1.96 * se), 4),
#   hi95   = round(pmin(1, resNewcomb$CPPP + 1.96 * se), 4),
#   row.names = NULL
# )


