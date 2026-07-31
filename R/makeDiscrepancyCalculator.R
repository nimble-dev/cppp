#' Build a discrepancy calculator
#'
#' Turns one or more discrepancy specifications and a simulation specification
#' into a function that computes discrepancy values from a model. For every
#' posterior draw it gives the discrepancies of the dataset you supply, and the
#' discrepancies of a replicate dataset simulated from that draw.
#'
#' The result is what [runCalibration()] and [runCalibrationNIMBLE()] take as
#' their `discFun` argument. The values are computed after the MCMC has run,
#' from its stored draws.
#'
#' Give the calculator its own copy of the model, `model$newModel()`, rather
#' than one an MCMC is using: it writes parameter values and data into the model
#' as it works, and puts them back when it is finished. The model must be
#' uncompiled — the calculator compiles its own copy.
#'
#' By default the discrepancies are compiled, which takes a little time once and
#' then runs far faster over the draws. Anything you write yourself has to be
#' code NIMBLE can compile; to use an R function, wrap it with
#' [nimble::nimbleRcall()]. Set `compile = FALSE` for a quick check without
#' waiting for compilation.
#'
#' Parameters may be derived nodes — `sigma`, say, when the model puts the prior
#' on `log(sigma)`. The drawn values are used as they are given.
#'
#' @param model A NIMBLE model.
#' @param discrepancies A [discrepancy()] specification, or a list of them.
#' @param simulation A [simulation()] specification.
#' @param paramNodes Character vector naming the model nodes to set from each
#'   posterior draw. These must appear among the column names of the draws.
#' @param compile Compile the model and the discrepancies? `TRUE` by default.
#'
#' @return A function of `(MCMCSamples, targetData, control, ...)`, returning a
#'   list of two matrices: `obs`, the discrepancies of `targetData`, and `sim`,
#'   those of the replicates. Both have one row per draw and one column per
#'   discrepancy, named after the discrepancies.
#'
#' @seealso [discrepancy()] and [simulation()] for writing the specifications,
#'   [runCalibrationNIMBLE()] for running a calibration.
#'
#' @examples
#' \dontrun{
#' calc <- makeDiscrepancyCalculator(
#'   model         = model$newModel(),
#'   discrepancies = list(discrepancy("mean"), discrepancy("deviance")),
#'   simulation    = simulation("conditional"),
#'   paramNodes    = c("mu", "sigma")
#' )
#'
#' d <- calc(MCMCSamples, targetData = y)
#' head(d$obs)
#' colMeans(d$sim >= d$obs)   # posterior predictive p-value per discrepancy
#' }
#' @export
makeDiscrepancyCalculator <- function(model, discrepancies, simulation, paramNodes,
                                      compile = TRUE) {

  ## A nimbleFunction is built against an uncompiled model, so we need one to
  ## start from even though the work then happens on the compiled copy.
  if (inherits(model, "CmodelBaseClass")) {
    stop("`model` must be an uncompiled NIMBLE model; the calculator compiles its own copy.",
         call. = FALSE)
  }

  discs   <- standardizeDiscrepancies(discrepancies)
  discs   <- lapply(discs, function(d) completeDiscrepancy(model, d))
  simSpec <- completeSimulation(model, simulation)

  ## SP: not sure if it is necessary to expand paramNodes here - it may be redundat
  paramNodes <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)

  discNames <- vapply(discs, function(d) d$name, character(1))
  if (anyDuplicated(discNames)) {
    stop("Each discrepancy needs its own name; repeated: ",
         paste(unique(discNames[duplicated(discNames)]), collapse = ", "),
         call. = FALSE)
  }

  ## SP: both discrepancy and simulation allows to specify dataNodes. For discrepancy
  ## those are the ones used in calculation and for simulation those nodes are the ones
  ## to simulate into. Nodes for discrepancy needs to be a subset (or match exactly)
  ## the ones in simulation
  for (d in discs) {
    unwritten <- setdiff(d$dataNodes, simSpec$dataNodes)
    if (length(unwritten) > 0L) {
      stop("Discrepancy '", d$name, "' reads data nodes the simulation does not set: ",
           paste(unwritten, collapse = ", "),
           ". Give `discrepancy()` and `simulation()` matching `dataNodes`.",
           call. = FALSE)
    }
  }

  ## SP: we get dependencies of parameters because we call calculate after changing values
  ## of the model parameters. self=FALSE to avoid overriding of deterministic nodes
  ## (e.g. user monitoring sigma (deterministic) instead of log_sigma (stochastic))
  paramDeps <- model$getDependencies(paramNodes, self = FALSE)
  ## After simulating a replicate: the data's own densities and anything below.
  dataDeps  <- model$getDependencies(simSpec$dataNodes, self = TRUE)
  restoreDeps <- model$topologicallySortNodes(unique(c(paramDeps, dataDeps)))

  discFuns <- lapply(discs, function(d) makeDiscrepancyNimbleFun(model, d))
  K <- length(discFuns)

  ## Compiled discrepancies read the compiled model's state, which is a
  ## different object from the R model. So once we compile, everything below
  ## works on the compiled copy: the values we write and the values the
  ## discrepancies read have to live in the same place.
  if (compile) {
    workModel <- compileNimble(model)
    discFuns  <- lapply(discFuns, function(f) compileNimble(f, project = model))
  } else {
    workModel <- model
  }

  function(MCMCSamples, targetData, control = NULL, ...) {

    MCMCSamples <- as.matrix(MCMCSamples)

    ## The columns need not be in the same order as `paramNodes`, so line them
    ## up by name. Matching against a matrix without names would give NA
    ## positions and read nonsense, so stop instead.
    absent <- setdiff(paramNodes, colnames(MCMCSamples))
    if (length(absent) > 0L) {
      stop("`MCMCSamples` has no column for: ", paste(absent, collapse = ", "),
           call. = FALSE)
    }
    paramCols <- match(paramNodes, colnames(MCMCSamples))

    if (!is.numeric(targetData) || length(targetData) != length(simSpec$dataNodes)) {
      stop("`targetData` must be a numeric vector with one value per data node (",
           length(simSpec$dataNodes), " expected).", call. = FALSE)
    }

    nDraws <- nrow(MCMCSamples)
    obs <- matrix(NA_real_, nDraws, K, dimnames = list(NULL, discNames))
    rep <- matrix(NA_real_, nDraws, K, dimnames = list(NULL, discNames))

    ## SP: on exit we need to leave the model as we found it.
    ## Restore once on exit is enough,
    ## because every draw overwrites these before reading them.
    savedData   <- values(workModel, simSpec$dataNodes)
    savedParams <- values(workModel, paramNodes)
    on.exit({
      values(workModel, simSpec$dataNodes) <- savedData
      values(workModel, paramNodes)        <- savedParams
      workModel$calculate(restoreDeps)
    })

    for (j in seq_len(nDraws)) {
      values(workModel, paramNodes)        <- MCMCSamples[j, paramCols]
      values(workModel, simSpec$dataNodes) <- targetData
      workModel$calculate(paramDeps)

      for (k in seq_len(K)) obs[j, k] <- discFuns[[k]]$run()

      workModel$simulate(simSpec$simulateNodes, includeData = TRUE)
      workModel$calculate(dataDeps)

      for (k in seq_len(K)) rep[j, k] <- discFuns[[k]]$run()
    }

    list(obs = obs, sim = rep)
  }
}
