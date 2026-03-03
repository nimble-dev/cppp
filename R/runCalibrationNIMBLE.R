#' Run calibration using a NIMBLE model
#'
#' @param model either an uncompiled or compiled nimbleModel with observed data set.
#' @param dataNames Optional character vector of data node names. If NULL,
#'   nodes flagged as data in the model are used.
#' @param paramNames Character vector of parameter node names to monitor. (SP: if NULL some default?)
#' @param MCMCSamples Optional matrix of posterior draws from the observed-data fit.
#'   If provided, `runCalibrationNIMBLE()` skips the main MCMC run and uses these
#'   draws as the posterior sample input to `runCalibration()`. The matrix must
#'   contain columns named by `paramNames`.
#' @param discFun Function `function(MCMCSamples, targetData, control)` that returns `list(obs, sim)` with one discrepancy value per posterior draw of `MCMCSamples`.
#' @param simulateNewDataFun Function `function(thetaRow, control)` that  simulates one replicated dataset from the posterior predictive. SP: We assume that new data is sampled from the posterior predictive of the model. In principle we may want to consider sampling from the prior predictive.
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
    paramNames,
    MCMCSamples  = NULL,
    discFun,
    simulateNewDataFun,
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
    message("Data nodes: ", paste(dataNames, collapse = ", "))
    message("Model class: ", paste(class(model), collapse = "/"))
    message("Compiled model class: ", paste(class(cmodel), collapse = "/"))
  }

  ## 1. Obtain posterior draws from observed data (either run MCMC or use user-supplied samples)
  if (is.null(MCMCSamples)) {

    ## Configure and compile MCMC for main chain
    if (is.null(mcmcConfFun)) {
      mcmcConfFun <- function(model) {
        configureMCMC(model, monitors = paramNames, print = FALSE)
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
  if (!all(paramNames %in% colnames(MCMCSamples))) {
    stop("paramNames missing from MCMCSamples: ",
         paste(setdiff(paramNames, colnames(MCMCSamples)), collapse = ", "))
  }

  ## Ensure we have a compiled MCMC object for replicated calibration runs
  if (!exists("cmcmc", inherits = FALSE)) {
    if (is.null(mcmcConfFun)) {
      mcmcConfFun <- function(model) {
        configureMCMC(model, monitors = paramNames, print = FALSE)
      }
    }
    mcmcConf       <- mcmcConfFun(model)
    mcmcUncompiled <- buildMCMC(mcmcConf)
    cmcmc          <- compileNimble(mcmcUncompiled, project = model, resetFunctions = TRUE)
  }

  ## Extract observed data from the model
  ## SP: need to think better if data is a vector/matrix/array - something else?
  observedData <- cmodel[[dataNames]]


  ## 4. Build MCMCFun for replicated datasets
  ## SP: do we want a function makeMCMCfun for specific configurations?
  MCMCFun <- function(targetData, control) {
    cmodel[[dataNames]] <- targetData
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
    disc = list(
      model      = model,
      dataNames  = dataNames,
      paramNames = paramNames
    ),
    draw = list()
  )

  control <- modifyList(defaultControl, control)

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

