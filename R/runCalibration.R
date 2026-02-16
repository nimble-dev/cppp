############################################################
## Generic calibration engine
############################################################

#' Run posterior predictive calibration
#'
#' Computes a calibrated posterior predictive p-value (CPPP) using
#' posterior predictive checks and replicated calibration worlds.
#'
#' The function is backend-agnostic: model fitting, data simulation,
#' and discrepancy evaluation are provided by user-supplied functions.
#'
#' @param MCMCSamples Matrix of posterior draws from the observed-data fit.
#'   Each row corresponds to one posterior draw of the model parameters.
#'
#' @param observedData Observed dataset. This dataset is treated
#'   as the conditioning data for the initial discrepancy calculation.
#'
#' @param MCMCFun Function defined as
#'   `function(targetData, control)` that runs an MCMC algorithm conditional on
#'   `targetData` and returns posterior samples as a matrix.
#'
#' @param simulateNewDataFun Function defined as
#'   `function(thetaRow, control)` that simulates one replicated dataset
#'   from the posterior predictive distribution given a single posterior
#'   draw `thetaRow`.
#'
#'
#' @param discFun Discrepancy extractor with signature
#'   `function(MCMCSamples, targetData, control)` returning a list with
#'   components `obs` and `sim`, each a numeric vector (one value per row of
#'   `MCMCSamples`).
#'
#'   In typical use, `discFun` is created by one of the package helpers:
#'   \describe{
#'     \item{\code{makeOfflineDiscFun()}}{Computes discrepancies by evaluating a
#'       user-specified discrepancy \eqn{D(data, \theta)} on `targetData` and on
#'       posterior predictive replicates simulated for each posterior draw.}
#'     \item{\code{makeColDiscFun()}}{Extracts precomputed discrepancy columns
#'       from `MCMCSamples` (online mode), returning them as `obs` and `sim`.}
#'   }
#'
#'   Custom `discFun` functions are also supported, as long as the return value
#'   has named components `obs` and `sim` of equal length.
#'
#' @param nReps Integer. Number of calibration replicates.
#'
#' @param drawIndexSelector Optional function
#'   `function(MCMCSamples, nReps, control)` returning a vector of row indices
#'   selecting which posterior draws are used as seeds for calibration worlds.
#'   If `NULL`, indices are chosen evenly over the posterior sample.
#'
#' @param control Optional list of arguments passed to components used by
#'   `runCalibration()`.
#'
#'   For convenience, `control` may be a flat list, in which case it is passed
#'   unchanged to all components.
#'
#'   Alternatively, `control` may be a structured list with named sublists:
#'   \describe{
#'     \item{mcmc}{Arguments passed to `MCMCFun`, typically controlling short
#'       MCMC runs in replicated calibration worlds (e.g., `niter`, `nburnin`).}
#'     \item{disc}{Arguments passed to `simulateNewDataFun` and `discFun`,
#'       typically including model objects, data node names, or parameter names.}
#'     \item{draw}{Arguments passed to `drawIndexSelector`, if provided.}
#'     \item{parallel}{Optional list controlling parallel execution of
#'       replicated calibration worlds. Supported fields:
#'       \describe{
#'         \item{workers}{Integer number of PSOCK workers. Default is 1 (serial).}
#'         \item{seed}{Optional integer seed for reproducible parallel RNG.}
#'         \item{export}{Optional character vector of object names (e.g. helper functions) to export to workers.}
#'         \item{packages}{Optional character vector of packages to load on workers.}
#'       }}
#'   }
#'
#' @param ... Not used currently.
#'
#' @return An object of class `cpppResult` containing:
#'   \describe{
#'     \item{CPPP}{Calibrated posterior predictive p-value.}
#'     \item{obsPPP}{Posterior predictive p-value for the observed data.}
#'     \item{repPPP}{Vector of posterior predictive p-values for each
#'       calibration world.}
#'     \item{discrepancies}{List containing observed and replicated discrepancies.}
#'     \item{drawnIndices}{Indices of posterior draws used as calibration seeds.}
#'   }
#'
#' @details
#' Calibration replicates are conditionally independent given the
#' posterior sample and may be executed in parallel using PSOCK clusters.
#' Objects referenced by `MCMCFun`, `discFun`, and `simulateNewDataFun` must
#' be serializable and available on worker processes.
#'
#' @export
#'

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

  ## 1. Setup and checks

  MCMCSamples <- as.matrix(MCMCSamples)
  nDraws      <- nrow(MCMCSamples)
  if (nDraws < 1L) stop("MCMCSamples must contain at least one row.")

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

  ##  Choose rows in MCMCSamples to simulate data for calibration
  if (is.null(drawIndexSelector)) {
    drawnIndices <- floor(seq(1, nDraws, length.out = nReps))
  } else {
    drawnIndices <- drawIndexSelector(MCMCSamples, nReps, drawControl)
  }
  if (length(drawnIndices) != nReps) {
    stop("drawIndexSelector must return exactly nReps indices.")
  }

  ## Set up parallelization
  parallelControl <- control$parallel
  workers <- if (is.list(parallelControl) && !is.null(parallelControl$workers)) {
    as.integer(parallelControl$workers)
  } else 1L

  seed <- if (is.list(parallelControl)) parallelControl$seed else NULL
  useParallel <- !is.na(workers) && workers > 1L

  ## Define all the rows of original MCMC (replication "seeds") that will be used to generate new replicates
  repSeeds <- MCMCSamples[drawnIndices, , drop = FALSE]

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

  ## function that run one replicate
  runOneRep <- function(r) {
    ## extract one row named vector from MCMC samples
    thetaRow <- repSeeds[r, , drop = TRUE]

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

    list(
      ppp  = mean(repSimDisc >= repObsDisc),
      disc = repDisc
    )

  }

  ## PARALLELIZATION
  if (!useParallel) {

    for (r in seq_len(nReps)) {
      out <- runOneRep(r)
      repPPP[r] <- out$ppp
      repDiscList[[r]] <- out$disc
    }

  } else {

    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    ## ---- PSOCK setup: load packages and export extra objects  ----
    parallelControl <- if (is.list(parallelControl)) parallelControl else list()

    pkgsToLoad <- parallelControl$packages
    if (!is.null(pkgsToLoad)) {
      pkgsToLoad <- as.character(pkgsToLoad)
      parallel::clusterEvalQ(cl, {
        for (p in pkgsToLoad) {
          if (!require(p, character.only = TRUE)) {
            stop(sprintf("Failed to load package '%s' on worker.", p))
          }
        }
        NULL
      })
    }

    extraExport <- parallelControl$export
    if (!is.null(extraExport)) {
      extraExport <- as.character(extraExport)
    } else {
      extraExport <- character(0)
    }

    ## export everything the workers need
    parallel::clusterExport(
      cl,
      varlist = unique(c(
        "repSeeds", "discControl", "mcmcControl",
        "simulateNewDataFun", "MCMCFun", "discFun", "runOneRep",
        extraExport
      )),
      envir = parent.frame()
    )

    ## To make runCalibrationNIMBLE work
    ## one-time worker initialization (optional)
    if (is.list(parallelControl) && is.function(parallelControl$init)) {
      parallel::clusterCall(cl, parallelControl$init)
    }

    ok <- unlist(parallel::clusterEvalQ(cl, exists("workerCtx", envir = .GlobalEnv)))
    if (!all(ok)) stop("PSOCK init ran, but workerCtx was not created on all workers.")

    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }

    outs <- parallel::parLapply(cl, X = seq_len(nReps), fun = runOneRep)

    repPPP <- vapply(outs, function(z) z$ppp, numeric(1))
    repDiscList <- lapply(outs, function(z) z$disc)
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




