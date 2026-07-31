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
#' as it works, and puts them back when it is finished.
#'
#' Parameters may be derived nodes — `sigma`, say, when the model puts the prior
#' on `log(sigma)`. The drawn values are used as they are given.
#'
#' @param model A NIMBLE model.
#' @param discrepancies A [discrepancy()] specification, or a list of them.
#' @param simulation A [simulation()] specification.
#' @param paramNodes Character vector naming the model nodes to set from each
#'   posterior draw. These must appear among the column names of the draws.
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
makeDiscrepancyCalculator <- function(model, discrepancies, simulation, paramNodes) {

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
    savedData   <- values(model, simSpec$dataNodes)
    savedParams <- values(model, paramNodes)
    on.exit({
      values(model, simSpec$dataNodes) <- savedData
      values(model, paramNodes)        <- savedParams
      model$calculate(restoreDeps)
    })

    for (j in seq_len(nDraws)) {
      values(model, paramNodes)        <- MCMCSamples[j, paramCols]
      values(model, simSpec$dataNodes) <- targetData
      model$calculate(paramDeps)

      for (k in seq_len(K)) obs[j, k] <- discFuns[[k]]$run()

      model$simulate(simSpec$simulateNodes, includeData = TRUE)
      model$calculate(dataDeps)

      for (k in seq_len(K)) rep[j, k] <- discFuns[[k]]$run()
    }

    list(obs = obs, sim = rep)
  }
}
