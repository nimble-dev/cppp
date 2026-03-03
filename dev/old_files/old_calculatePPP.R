###########
#' Compute the posterior predictive p
#'
#' function that simulate samples from the posterior predictive and calculate summaries
#'
#' @param data: original data
#' @param origMCMCSamples: matrix of posterior samples
#' @param ndraws: number of samples to simulate
#' @param sampleFromPP: function that can sample new data from the posterior predictive using the provided samples
#' @param discrepancy: a function that calculates the discrepancy (statistics) for the data (original or sampled)
#' @param ...: additional arguments to be passed
#'
#' @return list containing the p-value, sampled discrepancies and original one
#'
#' @export

calculatePPP <- function(data, origMCMCSamples, ndraws, sampleFromPP, discrepancy, ...) {
  ## Assumption:  one discrepancy with univariate output, do not depend on parameters

  if(ndraws > nrow(origMCMCSamples)) stop("want to simulate more samples than avaiable MCMC ones")
  nObs <- length(data)[1] ## TMP

  ppSamples <- list()
  discPP <- numeric(ndraws)

  ## simulate samples from the posterior predictive and calculate discrepancies
  for(i in 1:ndraws) {
    ppSamples[[i]] <- sampleFromPP(nObs = nObs, origMCMCSamplesRow = origMCMCSamples[i, ] )

    ## calculate the discrepancy for each sample
    discPP[i] <- discrepancy(ppSamples[[i]])

  }


  ## calculate the discrepancy for the original data
  discData <- discrepancy(data, ...)

  ## calculate the p-value
  ppp <- sum(discPP - discData >= 0) / ndraws

  ## return the p-value, sampled discrepancies and original one
  list(ppp = ppp, sampledDisc = discPP, discData = discData)

}
