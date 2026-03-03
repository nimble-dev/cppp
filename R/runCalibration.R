############################################################
## Generic calibration engine (backend-agnostic)
############################################################

#' Run posterior predictive calibration
#'
#' @param MCMCSamples Matrix of posterior draws from the observed-data fit.
#' @param observedData Observed dataset (any R object).
#' @param MCMCFun Function `function(targetData, control)` that runs a short MCMC conditional on `targetData`and returns posterior samples. `targetData` is the dataset treated as observed for the MCMC run.
#' @param simulateNewDataFun Function `function(thetaRow, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param discFun Function `function(MCMCSamples, targetData, control)` that returns a list with components `obs` and `sim`, each a numeric vector over posterior draws. `targetData` is the dataset treated as observed for the discrepancy calculation.
#' @param nReps Number of calibration replications.
#' @param drawIndexSelector Optional function `function(MCMCSamples, nReps, control)`
#'   returning the indices of rows to use as seeds for calibration.
#' @param control Optional list of specific arguments passed to the
#'   components used by `runCalibration()`. For convenience, `control` may be
#'   a flat list, in which case it is passed unchanged to all components.
#'
#'   Alternatively, `control` may be a structured list with named sublists,
#'   which are routed to different stages of the calibration:
#'   \describe{
#'     \item{mcmc}{Arguments or objects passed to `MCMCFun`, typically used
#'       to control short MCMC runs in replicated calibration worlds.}
#'     \item{disc}{Arguments or objects passed to `simulateNewDataFun` and
#'       `discFun`, typically used for data simulation and discrepancy
#'       calculation (e.g., model objects, node names).}
#'     \item{draw}{Arguments passed to `drawIndexSelector`, if provided.}
#'   }
#'
#'   This structure allows backend-specific state (such as NIMBLE models
#'   and compiled MCMC objects) to be cleanly separated while remaining
#'   backward compatible with simpler uses.#' @param ... Not used currently.
#'
#' @return a list (future `S3` class `cpppResult` objects) containing the CPPP, observed and replicated repPPP, observed discrepancies and replicated discrepancies.
#' @export

runCalibration <- function(
    MCMCSamples,
    observedData,
    MCMCFun,
    simulateNewDataFun,
    discFun,
    nReps,
    drawIndexSelector = NULL,
    control = list(),
    ...
) {

  ## Check that MCMCSamples is matrix
  if (!is.matrix(MCMCSamples)) {
    stop("MCMCSamples must be a numeric matrix with one row per posterior draw.")
  }

  if (!is.numeric(MCMCSamples)) {
    stop("MCMCSamples must be a numeric matrix.")
  }

  nDraws <- nrow(MCMCSamples)
  if (is.null(nDraws) || nDraws < 1L) {
    stop("MCMCSamples must contain at least one row.")
  }
  ## check what controls contains by role
  mcmcControl <- if (!is.null(control$mcmc)) control$mcmc else control
  discControl <- if (!is.null(control$disc)) control$disc else control
  drawControl <- if (!is.null(control$draw)) control$draw else control

  if (isTRUE(control$verbose)) {
    message("runCalibration(): routing")
    message("  names(control): ", paste(names(control), collapse = ", "))
    message("  names(mcmcControl): ", paste(names(mcmcControl), collapse = ", "))
    message("  names(discControl): ", paste(names(discControl), collapse = ", "))
    message("  names(drawControl): ", paste(names(drawControl), collapse = ", "))
  }

  ## check messages
  verbose <- isTRUE(control$verbose)
  ## for parallel
  # progressEvery <- control$progressEvery
  # if (!is.null(progressEvery)) progressEvery <- as.integer(progressEvery)

  ##  Choose rows in MCMCSamples to simulate data for calibration
  if (is.null(drawIndexSelector)) {
    drawnIndices <- floor(seq(1, nDraws, length.out = nReps))
  } else {
    drawnIndices <- drawIndexSelector(MCMCSamples, nReps, drawControl)
  }
  if (length(drawnIndices) != nReps) {
    stop("drawIndexSelector must return exactly nReps indices.")
  }

  ## 2. Discrepancies + PPP for the observed data
  obsDisc <- discFun(MCMCSamples  = MCMCSamples,
                     targetData   = observedData,
                     control      = discControl)
  if (!all(c("obs", "sim") %in% names(obsDisc))) {
    stop("discFun must return a list with components 'obs' and 'sim'.")
  }
  ## discrepancy computed on real data
  origObsDisc <- as.numeric(obsDisc$obs)
  origSimDisc <- as.numeric(obsDisc$sim)
  if (length(origObsDisc) != length(origSimDisc)) {
    stop("'obs' and 'sim' must have the same length.")
  }

  # scalar PPP for the observed data
  obsPPP <- mean(origSimDisc >= origObsDisc)

  ## 3. Calibration worlds: PPP for each replicates
  repPPP <- numeric(nReps)
  repDiscList <- vector("list", nReps)

  for (r in seq_len(nReps)) {
    ## extract one row named vector from MCMC samples
    thetaRow <- MCMCSamples[drawnIndices[r], , drop = TRUE]

    # 3a. simulate a new dataset y^(r) from posterior predictive of the original model
    newData <- simulateNewDataFun(thetaRow = thetaRow,
                                  control  = discControl)

    # 3b. fit model on y^(r) (short chain)
    repMCMC <- MCMCFun(targetData = newData,
                         control  = mcmcControl)

    # 3c. compute discrepancies + PPP in this world
    repDisc <- discFun(MCMCSamples = repMCMC,
                       targetData = newData,
                         control  = discControl)

    if (!all(c("obs", "sim") %in% names(repDisc))) {
      stop("discFun must return a list with components 'obs' and 'sim'.")
    }
    repObsDisc <- as.numeric(repDisc$obs)
    repSimDisc <- as.numeric(repDisc$sim)
    if (length(repObsDisc) != length(repSimDisc)) {
      stop("'obs' and 'sim' must have the same length in each calibration world.")
    }

    repPPP[r]         <- mean(repSimDisc >= repObsDisc)
    repDiscList[[r]] <- repDisc

  }

  ## 4. CPPP: how extreme obsPPP is under the calibration distribution
  CPPP <- mean(repPPP <= obsPPP)

  ## 5. Collect all discrepancies
  discrepancies <- list(
    obs = list(
      obs = origObsDisc,
      sim = origSimDisc
    ),
    rep = lapply(repDiscList, function(d) {
      list(
        obs = as.numeric(d$obs),
        sim = as.numeric(d$sim)
      )
    })
  )

  ## 6. Return cpppResult object
  newCpppResult(
    CPPP          = CPPP,
    repPPP        = repPPP,
    obsPPP        = obsPPP,
    discrepancies = discrepancies,
    drawnIndices  = drawnIndices
  )
}




