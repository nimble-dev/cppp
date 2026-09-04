#' Run calibration using a NIMBLE model
#'
#' @param model An uncompiled nimbleModel, initialized with observed data. The model, the MCMC and the discrepancy calculator are compiled together in one call.
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

  ## 0. The model, the MCMC and the discrepancy calculator are compiled
  ## together in one call further down, so we have to start from an uncompiled
  ## model: NIMBLE cannot add these to a project that already exists.
  if (inherits(model, "CmodelBaseClass")) {
    stop("`model` must be an uncompiled NIMBLE model. The model, the MCMC and ",
         "the discrepancy calculator are compiled together in one call.",
         call. = FALSE)
  }
  if (!inherits(model, "RmodelBaseClass")) {
    stop("Argument 'model' must be a nimbleModel.", call. = FALSE)
  }

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

  ## 1. Build the discrepancy calculator's pieces from specifications, if
  ## given. We compuile the calculator's nimbleFunction, the model and
  ## the MCMC once
  calcParts <- NULL
  simSpec   <- NULL
  if (!is.null(discrepancies)) {
    if (!is.null(discFun)) {
      stop("Give either `discrepancies` or `discFun`, not both.", call. = FALSE)
    }

    simSpec <- if (is.null(simulation)) simulation("conditional") else simulation

    ## The dataset written into the model, the one read back out, and the one
    ## the engine treats as observed all have to be the same nodes. Take
    ## `dataNodes` as the answer unless the specification says otherwise.
    if (is.null(simSpec$dataNodes)) simSpec$dataNodes <- dataNodes

    calcParts <- buildDiscrepancyCalculator(
      model         = model,
      discrepancies = discrepancies,
      simulation    = simSpec,
      paramNodes    = paramNodes
    )
  }

  if (is.null(calcParts) && is.null(discFun)) {
    stop("Supply `discrepancies` (recommended) or `discFun`.", call. = FALSE)
  }
  if (is.null(calcParts) && is.null(simulateNewDataFun)) {
    stop("Supply `simulateNewDataFun`, or `discrepancies` so it can be built for you.",
         call. = FALSE)
  }

  ## 2. Build the MCMC.
  if (is.null(mcmcConfFun)) {
    mcmcConfFun <- function(model) {
      configureMCMC(model, monitors = paramNodes, print = FALSE)
    }
  }
  mcmcConf       <- mcmcConfFun(model)
  mcmcUncompiled <- buildMCMC(mcmcConf)

  ## 3. One compile, for everything.
  toCompile <- list(model, mcmcUncompiled)
  if (!is.null(calcParts)) toCompile <- c(toCompile, list(calcParts$calcNF))

  compiled <- compileNimble(toCompile)
  cmodel   <- compiled[[1]]
  cmcmc    <- compiled[[2]]
  if (!is.null(calcParts)) calcParts$calcNF <- compiled[[3]]

  ## 4. Finish the two functions that need the compiled pieces.
  if (!is.null(calcParts)) {
    discFun <- wrapDiscrepancyCalculator(calcParts)

    if (is.null(simulateNewDataFun)) {
      ## Built on the compiled model, so simulating a replicate runs in
      ## compiled code. It also leaves the model holding the draw it simulated
      ## from, which is where that replicate's chain then starts.
      simulateNewDataFun <- makeSimulateNewDataFun(
        model      = cmodel,
        simulation = simSpec,
        paramNodes = paramNodes
      )
    }
  }

  if (verbose) {
    message("Data nodes: ", paste(dataNodes, collapse = ", "))
    message("Model class: ", paste(class(model), collapse = "/"))
    message("Compiled model class: ", paste(class(cmodel), collapse = "/"))
  }

  ## 5. Obtain posterior draws from observed data (either run MCMC or use user-supplied samples)
  if (is.null(MCMCSamples)) {

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

  ## SP: extract observed data as a numeric vector over expanded data nodes
  observedData <- values(cmodel, dataNodes)

  ## 6. Build MCMCFun for replicated datasets
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
