############################################################
## Newcomb example in pure R (no NIMBLE)
## Uses: runCalibration() + makeOfflineDiscFun() from cppp
############################################################
library(cppp)

## 1) Data
lightPath <- system.file("examples", "light.txt", package = "cppp")

yObs <- read.table(lightPath)$V1
n <- length(yObs)

dataNames  <- "y"
paramNames <- c("mu", "sigma")

## 2) Posterior sampler for Normal model with prior p(mu, sigma) ∝ 1/sigma
## Returns matrix with columns mu, sigma (named), one row per draw.
samplePosteriorNormalJeffreys <- function(y, nDraws) {
  ybar <- mean(y)
  s2   <- var(y)                 # sample variance with denominator (n-1)
  df   <- length(y) - 1

  ## sigma^2 | y ~ Inv-chi-square(df, s2)
  sigma2 <- (df * s2) / stats::rchisq(nDraws, df = df)
  sigma  <- sqrt(sigma2)

  ## mu | sigma, y ~ Normal(ybar, sigma/sqrt(n))
  mu <- stats::rnorm(nDraws, mean = ybar, sd = sigma / sqrt(length(y)))

  out <- cbind(mu = mu, sigma = sigma)
  out
}

## 3) MCMCFun: "short chain" on targetData (here: direct posterior sampling)
## Signature matches runCalibration(): function(targetData, control)
MCMCFun <- function(targetData, control) {
  nDraws <- control$nDraws
  samplePosteriorNormalJeffreys(targetData, nDraws = nDraws)
}

## 4) simulateNewDataFun: posterior predictive simulator given one draw thetaRow
## Signature: function(thetaRow, control) where thetaRow is named vector
simulateNewDataFun <- function(thetaRow, control) {
  n <- control$n
  mu <- unname(thetaRow["mu"])
  sigma <- unname(thetaRow["sigma"])
  rnorm(n, mean = mu, sd = sigma)
}

## 5) Discrepancy (same structure as your asymmetry example)
asymmDisc <- function(data, thetaRow, control) {
  mu <- unname(thetaRow["mu"])
  dataSorted <- sort(data)
  abs(dataSorted[61] - mu) - abs(dataSorted[6] - mu)
}

## 6) Build offline discFun using package helper
discConfig <- list(
  simulateNewDataFun = simulateNewDataFun,
  discrepancy        = asymmDisc
)
discFun <- makeOfflineDiscFun(discConfig)

## 7) Run generic calibration
## Use structured control so routing is exercised:
control <- list(
  verbose = TRUE,
  mcmc = list(nDraws = 100),    # size of each posterior sample used in discFun + refits
  disc = list(n = n)            # used by simulateNewDataFun + discrepancy (if needed)
)

set.seed(1)

## main posterior draws from observed data
MCMCSamples <- samplePosteriorNormalJeffreys(yObs, nDraws = control$mcmc$nDraws)

res <- runCalibration(
  MCMCSamples        = MCMCSamples,
  observedData       = yObs,
  MCMCFun            = MCMCFun,
  simulateNewDataFun = simulateNewDataFun,
  discFun            = discFun,
  nReps              = 10,
  drawIndexSelector  = NULL,
  control            = control
)

print(res$CPPP)
print(res$repPPP)
print(res$obsPPP)

## Optional quick look at observed discrepancies
## obsDisc <- res$discrepancies$obs
## plot(obsDisc$obs, obsDisc$sim); abline(0,1)
