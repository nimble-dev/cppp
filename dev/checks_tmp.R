library(cpppPrototype)
####################
## Silly testing ##
####################
## matrix of posterior samples
## 2 samples for 2 variables
samples <- matrix(c(rnorm(100), rgamma(100, 1, 2)), 100,2)
colnames(samples) <- c("mu", "s2")

## function that can sample from the posterior predictive
## number of samples
## other arguments
## output is a list containing each posterior predictive sample
sampleFromPP <- function(nObs, origMCMCSamplesRow, ... ) {
  rnorm(nObs, mean = origMCMCSamplesRow["mu"], sd = sqrt(origMCMCSamplesRow["s2"]))
}

## function that calculates the discrepancy
## input is a vector of data
## output is a scalar
discrepancy <- function(data) {
  mean(data)
}

## test

data <- rnorm(10)
res <- calculatePPP(data, origMCMCSamples = samples, ndraws = 100, sampleFromPP, discrepancy)
hist(res$sampledDisc, breaks= 20)
abline(v = res$discData, col = "red")


##########################
sampleFromPosterior <- function(nSamples, data, params, ... ) {
  res <- matrix(NA, nSamples, length(params))
  colnames(res) <- params
  ## some initial starting value
  res[1, ] <- c(2, 1)
  for(i in 2:nSamples){
    res[i, params[1]]<- rnorm(1, mean = res[i-1, params[1]] + mean(data), sd = 1)
    res[i, params[2]]<- rgamma(1, 1, 2)
  }

  res
}


# nCalibrationReplicates = 10
# origMCMCSamples = samples
# newData <- sampleFromPP(nObs = 10, origMCMCSamples[1, ])
# sampleFromPosterior(nSamples = 10, data, params = c("mu", "s2"))


calibratedSamples <- runCalibration(data = data,
               nCalibrationReplicates = 100,
               MCMCcontrol = list(niter = 200,
                                  thin = 1,
                                  nburnin = 0),
               sampleFromPosterior,
               sampleFromPP,
               discrepancy,
               origMCMCSamples = samples)

hist(calibratedSamples, breaks= 20)


calculateCPPP(data = data,
              nCalibrationReplicates = 100,
              MCMCcontrol = list(niter = 200,
                                 thin = 1,
                                 nburnin = 0),
              sampleFromPosterior,
              sampleFromPP,
              discrepancy,
              origMCMCSamples = samples)
