############################################################
## Generic calibration engine (backend-agnostic)
############################################################

#' Run posterior predictive calibration
#'
#' @param MCMCSamples Matrix of posterior draws from the observed-data fit.
#' @param observedData Observed dataset (any R object).
#' @param MCMCFun Function `function(targetData, control)` that runs a short MCMC conditional on `targetData`and returns posterior samples. `targetData` is the dataset treated as observed for the MCMC run.
#' @param simulateNewDataFun Function `function(thetaRow, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param discFun Function `function(MCMCSamples, targetData, control)` that
#'   returns a list with components `obs` and `sim`. `targetData` is the dataset
#'   treated as observed for the discrepancy calculation.
#'
#'   Design note: `discFun` is the one place where the discrepancy values are
#'   produced. There are two ways to produce them, and `discFun` hides which one
#'   is used:
#'   \itemize{
#'     \item offline: it computes the discrepancies here in R, for each draw,
#'       on `targetData` and on a simulated replicate `y*`.
#'     \item online: the discrepancies were already computed during the MCMC, so
#'       it just reads them back from columns of `MCMCSamples`.
#'   }
#'   Either way, `runCalibration()` receives the same thing and does not need to
#'   know how the numbers were made. Users are not expected to write `discFun`
#'   by hand: it will be built automatically from a discrepancy/simulation
#'   specification. The specification is what users work with; `discFun` is an
#'   internal detail.
#'
#'   `obs` and `sim` each hold one discrepancy value per posterior draw. With a
#'   single discrepancy this is a plain numeric vector (one value per draw).
#'   With several discrepancies it is an `nDraws x K` matrix: one row per draw,
#'   one column per discrepancy. `obs` and `sim` must have the same shape.
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
  ## Normalize the discFun output to nDraws x K matrices (one column per
  ## discrepancy). The discrepancy itself is always a scalar per (draw,
  ## discrepancy); the matrix is just those scalars stacked. We accept a plain
  ## vector (K = 1) for backward compatibility and so the single- and
  ## multiple-discrepancy cases share one downstream code path. See the
  ## `discFun` design note in this function's docs.
  origObs <- asDiscMatrix(obsDisc$obs, obsDisc)
  origSim <- asDiscMatrix(obsDisc$sim, obsDisc)
  checkDiscMatrices(origObs, origSim)

  discNames <- colnames(origObs)
  K <- ncol(origObs)

  # per-discrepancy PPP for the observed data (named length-K vector)
  obsPPP <- colMeans(origSim >= origObs)
  names(obsPPP) <- discNames

  ## 3. Calibration worlds: per-discrepancy PPP for each replicate
  repPPP <- matrix(NA_real_, nrow = nReps, ncol = K,
                   dimnames = list(NULL, discNames))
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

    repObs <- asDiscMatrix(repDisc$obs, repDisc, discNames)
    repSim <- asDiscMatrix(repDisc$sim, repDisc, discNames)
    checkDiscMatrices(repObs, repSim, K)

    repPPP[r, ]      <- colMeans(repSim >= repObs)
    repDiscList[[r]] <- list(obs = repObs, sim = repSim)

  }

  ## 4. CPPP: how extreme obsPPP is under the calibration distribution,
  ## computed separately for each discrepancy (named length-K vector).
  CPPP <- vapply(seq_len(K),
                 function(k) mean(repPPP[, k] <= obsPPP[k]),
                 numeric(1))
  names(CPPP) <- discNames

  ## 5. Collect all discrepancies
  discrepancies <- list(
    names = discNames,
    obs = list(
      obs = origObs,
      sim = origSim
    ),
    rep = repDiscList
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


#' Coerce discrepancy output to an nDraws x K matrix
#'
#' A `discFun` may return one discrepancy per draw (a numeric vector) or several
#' (a matrix with one column per discrepancy). This helper normalizes both to a
#' matrix with named columns so the engine can treat the single- and
#' multiple-discrepancy cases uniformly.
#'
#' @param x Numeric vector or matrix returned by a `discFun` component.
#' @param disc The full list returned by the `discFun`, inspected for a `names`
#'   element giving discrepancy names.
#' @param fallbackNames Optional character vector of column names to use when
#'   none are available from `x` or `disc`.
#'
#' @return A numeric matrix with named columns.
#' @keywords internal
asDiscMatrix <- function(x, disc = NULL, fallbackNames = NULL) {
  m <- if (is.matrix(x)) {
    storage.mode(x) <- "double"
    x
  } else {
    matrix(as.numeric(x), ncol = 1L)
  }

  if (is.null(colnames(m))) {
    nm <- fallbackNames
    if (is.null(nm) && !is.null(disc$names)) nm <- disc$names
    if (is.null(nm)) {
      nm <- if (ncol(m) == 1L) "discrepancy" else paste0("discrepancy", seq_len(ncol(m)))
    }
    colnames(m) <- nm
  }

  m
}


#' Validate paired discrepancy matrices
#'
#' @param obs Numeric matrix of observed-side discrepancies.
#' @param sim Numeric matrix of replicated-side discrepancies.
#' @param K Optional expected number of discrepancies (columns).
#' @keywords internal
checkDiscMatrices <- function(obs, sim, K = NULL) {
  if (!all(dim(obs) == dim(sim))) {
    stop("discFun 'obs' and 'sim' must have the same dimensions.", call. = FALSE)
  }
  if (!is.null(K) && ncol(obs) != K) {
    stop("discFun returned a different number of discrepancies across worlds.",
         call. = FALSE)
  }
  invisible(TRUE)
}




