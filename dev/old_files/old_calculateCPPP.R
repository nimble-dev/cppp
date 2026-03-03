####################################################
#' Run calibration
#'
#' This function runs the calibration (aka parametric bootstrap) to obtain the calibrated posterior predictive p-value.
#'
#' @param data: original data
#' @param nCalibrationReplicates: number of calibration replicates
#' @param MCMCcontrol: list with number of posterior samples per replicates (and optional thinning and burnin)
#' @param sampleFromPosterior: function that can sample from the posterior distribution
#' @param sampleFromPP: function that can sample new data from the posterior predictive using some provided samples
#' @param discrepancy: a function that calculates the discrepancy (statistics) for the data (original or sampled)
#' @param origMCMCSamples: matrix of posterior samples from the original data
#' @param ...: additional arguments to be passed
#'
#'
#' @export

runCalibration <- function(data,
                           nCalibrationReplicates = 100,
                           MCMCcontrol = list(niter = 200,
                                              thin = 1,
                                              nburnin = 0),
                           sampleFromPosterior,
                           sampleFromPP,
                           discrepancy,
                           origMCMCSamples,
                           ...) {
  ## run the calibration:
  ## 1. sample new data from the posterior predictive
  ## 2. treat the new data as the observed data and get the posterior distribution
  ## 3. repeat the above steps nCalibrationReplicates times

  paramsNames = colnames(origMCMCSamples)

  ## output - vector of PPPs
  calibration <- numeric(nCalibrationReplicates)

  ## choose MCMC iteration indices evenly spaces (SP: alternatively at random)
  rowsToUse <- floor( seq(1, nrow(origMCMCSamples), length=nCalibrationReplicates) )

  ## loop over the replicates
  for(i in 1:nCalibrationReplicates) {
    ## sample new data from the posterior predictive
    ## using random mcmc iteration
    newData <- sampleFromPP(1, origMCMCSamples[rowsToUse[i], ])

    ## treat the new data as the observed data and get the posterior distribution
    replicateMCMCSamples <- sampleFromPosterior(nSamples = MCMCcontrol$niter,
                                                data = newData,
                                                params = paramsNames)

    ## Calculate PPP
    resCalibrate <- calculatePPP(data = newData,
                                 origMCMCSamples = replicateMCMCSamples,
                                 ndraws = MCMCcontrol$niter,
                                 sampleFromPP, discrepancy)

    calibration[i] <- resCalibrate$ppp
  }

  ## return the results
  calibration
}

#######################################
#' Calculate the calibrated posterior predictive p-values
#'
#' function that calculate the cppp samples from the posterior predictive
#' and calculate the cppp
#'
#' @param data: original data
#' @param nCalibrationReplicates: number of calibration replicates
#' @param MCMCcontrol: list with number of posterior samples per replicates (and optional thinning and burnin)
#' @param sampleFromPosterior: function that can sample from the posterior distribution
#' @param sampleFromPP: function that can sample new data from the posterior predictive using some provided samples
#' @param discrepancy: a function that calculates the discrepancy (statistics) for the data (original or sampled)
#' @param origMCMCSamples: matrix of posterior samples from the original data
#' @param ...: additional arguments to be passed
#'
#'
#' @export

calculateCPPP <- function(data,
                          nCalibrationReplicates = 100,
                          MCMCcontrol = list(niter = 200,
                                             thin = 1,
                                             nburnin = 0),
                          sampleFromPosterior,
                          sampleFromPP,
                          discrepancy,
                          origMCMCSamples,
                          ...) {

  #### Resume here
  ## calculatePPP on the original data if not provided
  pppOriginal <- calculatePPP(data, origMCMCSamples,
                              ndraws = dim(origMCMCSamples)[1],
                              sampleFromPP,
                              discrepancy)

  ## run the calibration
  pppSamples <- runCalibration(data = data,
                               nCalibrationReplicates = 100,
                               MCMCcontrol = list(niter = 200,
                                                  thin = 1,
                                                  nburnin = 0),
                               sampleFromPosterior,
                               sampleFromPP,
                               discrepancy,
                               origMCMCSamples)

  ## calculate the cppp
  cppp <- sum(pppSamples - pppOriginal$ppp <= 0) / nCalibrationReplicates

  ## return the cppp
  cppp
}
