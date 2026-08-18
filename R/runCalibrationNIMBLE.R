#' Run calibration using a NIMBLE model
#'
#' @param model either an uncompiled or compiled nimbleModel, initialized with observed data.
#' @param dataNames Optional character vector of data node names. If NULL,
#'   nodes flagged as data in the model are used.
#' @param paramNames Optional character vector of parameter node names to monitor.
#' If NULL, all stochastic non-data nodes in the model are used.

#' @param MCMCSamples Optional matrix of posterior draws from the observed-data fit.
#'   If provided, `runCalibrationNIMBLE()` skips the main MCMC run and uses these
#'   draws as the posterior sample input to `runCalibration()`. The matrix must contain columns
#'   corresponding to the expanded parameter nodes from `paramNames`.
#' @param discrepancies A [discrepancy()] specification, or a list of them. Give
#'   these and the discrepancy calculator is built for you, instead of passing
#'   `discFun`.
#' @param simulation A [simulation()] specification, saying how a replicate
#'   dataset is generated. Defaults to `simulation("conditional")`.
#' @param discFun Function `function(MCMCSamples, targetData, control)` that returns `list(obs, sim)` with one discrepancy value per posterior draw of `MCMCSamples`. Supply this or `discrepancies`, not both.
#' @param simulateNewDataFun Function `function(thetaRow, control)` that  simulates one replicated dataset from the posterior predictive. Optional when `discrepancies` is given: it is then built from `simulation`. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
#' @param nReps Number of calibration replications.
#' @param MCMCcontrolMain List with `niter`, `nburnin`, `thin` for main chain.
#' @param MCMCcontrolRep List with `niter`, `nburnin`, `thin` for calibration chains.
#' @param mcmcConfFun Optional function `function(model)` returning an MCMC configuration.
#' @param drawIndexSelector Optional row selector (see runCalibration()).
#' @param control List of additional options passed to discrepancy / helpers.
#' @param ... Not used currently.
#' @export

runCalibrationNIMBLE <- function(
    model,
    dataNames    = NULL,
    paramNames   = NULL,
    MCMCSamples  = NULL,
    discrepancies = NULL,
    simulation    = NULL,
    discFun = NULL,
    simulateNewDataFun = NULL,
    nReps       = 100,
    MCMCcontrolMain = list(niter = 5000, nburnin = 1000, thin = 1),
    MCMCcontrolRep  = list(niter = 500,  nburnin = 0,    thin = 1),
    mcmcConfFun = NULL,
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
  dataNodes <- model$expandNodeNames(dataNames, returnScalarComponents = TRUE)
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
  paramNodes <- model$expandNodeNames(paramNames, returnScalarComponents = TRUE)
  if (length(paramNodes) == 0) {
    stop("paramNames did not match to any stochastic non-data nodes.")
  }

  ## 0. Build the discrepancy calculator from specifications, if given.
  if (!is.null(discrepancies)) {
    if (!is.null(discFun)) {
      stop("Give either `discrepancies` or `discFun`, not both.", call. = FALSE)
    }
    ## Named apart from the `simulation` argument so it is clear which is which.
    simSpec <- if (is.null(simulation)) simulation("conditional") else simulation

    ## The dataset written into the model, the one read back out, and the one
    ## the engine treats as observed all have to be the same nodes. Take
    ## `dataNames` as the answer unless the specification says otherwise.
    if (is.null(simSpec$dataNodes)) simSpec$dataNodes <- dataNames

    ## The calculator needs its own copy of the model
    baseModel <- if (inherits(model, "CmodelBaseClass")) model$Rmodel else model

    discFun <- makeDiscrepancyCalculator(
      model         = baseModel$newModel(),
      discrepancies = discrepancies,
      simulation    = simSpec,
      paramNodes    = paramNodes
    )

    ## And the replicate-simulating function, from the same specification.
    if (missing(simulateNewDataFun) || is.null(simulateNewDataFun)) {
      simulateNewDataFun <- makeSimulateNewDataFun(
        model      = baseModel$newModel(),
        simulation = simSpec,
        paramNodes = paramNodes
      )
    }
  }

  if (is.null(discFun)) {
    stop("Supply `discrepancies` (recommended) or `discFun`.", call. = FALSE)
  }
  if (is.null(simulateNewDataFun)) {
    stop("Supply `simulateNewDataFun`, or `discrepancies` so it can be built for you.",
         call. = FALSE)
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

  if (verbose) {
    message("Data nodes: ", paste(dataNodes, collapse = ", "))
    message("Model class: ", paste(class(model), collapse = "/"))
    message("Compiled model class: ", paste(class(cmodel), collapse = "/"))
  }

  ## 1. Obtain posterior draws from observed data (either run MCMC or use user-supplied samples)
  if (is.null(MCMCSamples)) {

    ## Configure and compile MCMC for main chain
    if (is.null(mcmcConfFun)) {
      mcmcConfFun <- function(model) {
        configureMCMC(model, monitors = paramNodes, print = FALSE)
      }
    }
    mcmcConf       <- mcmcConfFun(model)
    mcmcUncompiled <- buildMCMC(mcmcConf)
    cmcmc          <- compileNimble(mcmcUncompiled, project = model, resetFunctions = TRUE)

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
    mcmcConf       <- mcmcConfFun(model)
    mcmcUncompiled <- buildMCMC(mcmcConf)
    cmcmc          <- compileNimble(mcmcUncompiled, project = model, resetFunctions = TRUE)
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
    message("modifint the control for runCalibration:")
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
      model      = model,
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

  ## 5. call the generic runCalibration function
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
