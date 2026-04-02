#' Run calibration using a NIMBLE model
#'
#' @param model either an uncompiled or compiled nimbleModel with observed data set.
#' @param dataNames Optional character vector of data node names. If NULL,
#'   nodes flagged as data in the model are used.
#' @param paramNames Optional character vector of parameter node names to monitor.
#' If NULL, all stochastic non-data nodes in the model are used.

#' @param MCMCSamples Optional matrix of posterior draws from the observed-data fit.
#'   Each row corresponds to one posterior draw of the model parameters and columns names should contains `paramNames`.
#'   If provided, `runCalibrationNIMBLE()` skips the main MCMC run and uses these
#'   draws as the posterior sample input to `runCalibration()`. The matrix must contain columns
#'   corresponding to the expanded parameter nodes from `paramNames`.
#' @param discFun Function `function(MCMCSamples, targetData, control)` that returns `list(obs, sim)` with one discrepancy value per posterior draw of `MCMCSamples`.
#' @param simulateNewDataFun Function `function(thetaRow, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param MCMCcontrolMain List with `niter`, `nburnin`, `thin` for main chain.
#' @param MCMCcontrolRep List with `niter`, `nburnin`, `thin` for calibration chains.
#' @param mcmcConfFun Optional function `function(model)` returning an MCMC configuration.
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
#'   }#' @param ... Not used currently.
#' @export

runCalibrationNIMBLE <- function(
    model,
    dataNames    = NULL,
    paramNames   = NULL,
    MCMCSamples  = NULL,
    discFun,
    simulateNewDataFun,
    MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1),
    MCMCcontrolRep  = list(niter = 500,  nburnin = 0,    thin = 1),
    mcmcConfFun = NULL,
    nReps       = 100,
    drawIndexSelector = NULL,
    control = list(),
    ...
) {

  verbose <- isTRUE(control$verbose)

  ## 0. Data names and checks
  ## if dataNames is not provided, then use all nodes in the model that are data
  if (is.null(dataNames)) {
    dataNames <- model$getNodeNames(dataOnly = TRUE)
  }
  ## expand to nodes
  dataNodes <- model$expandNodeNames(dataNames)
  # ensure dataNames correspond to stochastic nodes
  testDataNames <- all(dataNodes %in%
                         model$getNodeNames(stochOnly = TRUE))
  if (!testDataNames) {
    stop("All dataNames must be stochastic nodes in the model.")
  }

  ## 0. deal with paramNodes. If missing we take all the stochastic nodes that are not data
  if (is.null(paramNames)) {
    paramNames <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  }
  paramNodes <- model$expandNodeNames(paramNames)
  if (length(paramNodes) == 0) {
    stop("paramNames did not match to any stochastic non-data nodes.")
  }

  ## check if the model is compiled model
  if (inherits(model, "CmodelBaseClass")) {
    ## Model is already compiled
    cmodel <- model
  } else if (inherits(model, "RmodelBaseClass")) {
    ## Model is not compiled yet
    cmodel <- compileNimble(model)
  } else {
    stop("Argument 'model' must be a nimbleModel or a compiled nimble model.")
  }

  ## ensure we always have an R model for worker initialization
  rModel <- if (inherits(model, "CmodelBaseClass")) model$getModel() else model

  if (verbose) {
    message("Data nodes: ", paste(dataNodes, collapse = ", "))
    message("Model class: ", paste(class(model), collapse = "/"))
    message("Compiled model class: ", paste(class(cmodel), collapse = "/"))
  }

  ## -------------------------------------------------------------------------------------------- ##
  ## 1. Obtain posterior draws from observed data (either run MCMC or use user-supplied samples)
  ## -------------------------------------------------------------------------------------------- ##
  if (is.null(MCMCSamples)) {

    ## Configure and compile MCMC for main chain
    if (is.null(mcmcConfFun)) {
      mcmcConfFun <- function(model) {
        configureMCMC(model, monitors = paramNodes, print = FALSE)
      }
    }
    mcmcConf       <- mcmcConfFun(rModel)
    mcmcUncompiled <- buildMCMC(mcmcConf)
    cmcmc          <- compileNimble(mcmcUncompiled, project = rModel, resetFunctions = TRUE)

    ## Run main chain on observed data
    obsMCMC <- runMCMC(
      cmcmc,
      niter   = MCMCcontrolMain$niter,
      nburnin = MCMCcontrolMain$nburnin,
      thin    = MCMCcontrolMain$thin
    )
    MCMCSamples <- as.matrix(obsMCMC)

    if (verbose) {
      message("Main MCMC finished")
      message("MCMCSamples dim: ", paste(dim(MCMCSamples), collapse = " x "))
      message("MCMCSamples columns: ", paste(colnames(MCMCSamples), collapse = ", "))
      message("paramNames: ", paste(paramNames, collapse = ", "))
    }

  } else {

    ## User supplied posterior draws
    MCMCSamples <- as.matrix(MCMCSamples)

    if (verbose) {
      message("Using user-supplied MCMCSamples")
      message("MCMCSamples dim: ", paste(dim(MCMCSamples), collapse = " x "))
      message("MCMCSamples columns: ", paste(colnames(MCMCSamples), collapse = ", "))
      message("paramNames: ", paste(paramNames, collapse = ", "))
    }
  }

  ## Validate samples contain required params
  if (!all(paramNodes %in% colnames(MCMCSamples))) {
    stop("paramNames missing from MCMCSamples: ",
         paste(setdiff(paramNodes, colnames(MCMCSamples)), collapse = ", "))
  }

    ## Ensure we have a compiled MCMC object for replicated calibration runs
  if (!exists("cmcmc", inherits = FALSE)) {
    if (is.null(mcmcConfFun)) {
      mcmcConfFun <- function(model) {
        configureMCMC(model, monitors = paramNodes, print = FALSE)
      }
    }
    mcmcConf       <- mcmcConfFun(rModel)
    mcmcUncompiled <- buildMCMC(mcmcConf)
    cmcmc          <- compileNimble(mcmcUncompiled, project = rModel, resetFunctions = TRUE)
  }

  ## SP: extract observed data as a numeric vector over expanded data nodes
  observedData <- values(cmodel, dataNodes)

  ## 4. Build MCMCFun for replicated datasets
  MCMCFun <- function(targetData, control) {
    if (!is.numeric(targetData) || length(targetData) != length(dataNodes)) {
      stop("targetData must be a numeric vector with one entry per data node.")
    }

    values(cmodel, dataNodes) <- targetData

    repMCMC <- runMCMC(
      cmcmc,
      niter   = control$niter,
      nburnin = control$nburnin,
      thin    = control$thin
    )
    as.matrix(repMCMC)
  }

  if (verbose) {
    message("modifiying the control for runCalibration:")
    message("  mcmc fields: ", paste(names(control$mcmc), collapse = ", "))
    message("  disc fields: ", paste(names(control$disc), collapse = ", "))
    message("  draw fields: ", paste(names(control$draw), collapse = ", "))
  }

  defaultControl <- list(
    mcmc = MCMCcontrolRep,
    ## SP: it may be helful to pass the nodes
    ## but depends on defaults for the discrepancy
    ## functions
    disc = list(
      model      = rModel,
      dataNames  = dataNames,
      dataNodes  = dataNodes,
      paramNames = paramNames,
      paramNodes = paramNodes
    ),
    # disc = list(
    #   model      = model,
    #   dataNames  = dataNames,
    #   paramNames = paramNames
    # ),
    draw = list()
  )

  control <- modifyList(defaultControl, control)

  ## -------------------------------------------------------------------------------------------- ##
  ## SP - NEW at 20260202: Set up for parallelization
  ## -------------------------------------------------------------------------------------------- ##
  workers <- if (is.list(control$parallel) && !is.null(control$parallel$workers)) {
    as.integer(control$parallel$workers)
  } else 1L
  useParallel <- !is.na(workers) && workers > 1L

  if (useParallel) {
    ## ensure nimble is loaded on workers
    if (is.null(control$parallel$packages)) control$parallel$packages <- character(0)
    control$parallel$packages <- unique(c(control$parallel$packages, "nimble"))

    ## ensure needed objects/functions are exported to workers
    if (is.null(control$parallel$export)) control$parallel$export <- character(0)
    control$parallel$export <- unique(c(
      control$parallel$export,
      "rModel", "dataNames", "paramNames", "mcmcConfFun"
    ))

    ## function to initialize each worker: make a worker-local compiled context
    control$parallel$init <- function() {
      tryCatch({
        workerCtx <- new.env(parent = emptyenv())
        workerCtx$repModel  <- rModel$newModel(replicate = TRUE)
        workerCtx$cRepModel <- nimble::compileNimble(workerCtx$repModel)

        conf <- if (is.null(mcmcConfFun)) {
          nimble::configureMCMC(workerCtx$repModel, monitors = paramNames, print = FALSE)
        } else {
          mcmcConfFun(workerCtx$repModel)
        }
        mcmcUncompiled <- nimble::buildMCMC(conf)
        workerCtx$cRepMCMC <- nimble::compileNimble(mcmcUncompiled, project = workerCtx$repModel, resetFunctions = TRUE)

        assign("workerCtx", workerCtx, envir = .GlobalEnv)
        TRUE
      }, error = function(e) {
        assign("workerCtx_init_error", conditionMessage(e), envir = .GlobalEnv)
        FALSE
      })
    }

  }

  ##### MCMFun

  if (!useParallel) {

    MCMCFun <- function(targetData, control) {
      cmodel[[dataNames]] <- targetData
      as.matrix(nimble::runMCMC(
        cmcmc,
        niter   = control$niter,
        nburnin = control$nburnin,
        thin    = control$thin
      ))
    }

  } else {

    MCMCFun <- function(targetData, control) {
      ctx <- get("workerCtx", envir = .GlobalEnv)
      ctx$cRepModel[[dataNames]] <- targetData
      as.matrix(nimble::runMCMC(
        ctx$cRepMCMC,
        niter   = control$niter,
        nburnin = control$nburnin,
        thin    = control$thin
      ))
    }

  }


  if (useParallel) {
    userSim  <- simulateNewDataFun
    userDisc <- discFun

    simulateNewDataFun <- function(thetaRow, control) {
      ctx <- get("workerCtx", envir = .GlobalEnv)
      control$model <- ctx$cRepModel   ## or ctx$repModel if user wants R model methods
      userSim(thetaRow, control)
    }

    discFun <- function(MCMCSamples, targetData, control) {
      ctx <- get("workerCtx", envir = .GlobalEnv)
      control$model <- ctx$cRepModel
      userDisc(MCMCSamples, targetData, control)
    }
  }


  ## 5. Call generic engine
  if (verbose) {
    message("Calling runCalibration() with nReps = ", nReps)
  }

  runCalibration(
    MCMCSamples  = MCMCSamples,
    observedData = observedData,
    MCMCFun      = MCMCFun,
    simulateNewDataFun  = simulateNewDataFun,
    discFun      = discFun,
    nReps        = nReps,
    drawIndexSelector  = drawIndexSelector,
    control       = control
  )

}

