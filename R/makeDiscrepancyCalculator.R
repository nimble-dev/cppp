#' Build a discrepancy calculator
#'
#' Combine what the user wrote — one or more [discrepancy()] specifications and
#' a [simulation()] specification — into a single function that computes
#' discrepancy values from the model, for the real data and for a replicate,
#' one pair per posterior draw.
#'
#' This is the calculating route: the values are worked out here, after the
#' MCMC has run. The other route is to have the MCMC compute them as it goes
#' and read the columns back afterwards. Either way the values are handed on to
#' [runCalibration()], which is what turns them into a PPP.
#'
#' Everything that can be prepared once is prepared here: the node names are
#' filled in from the model, and one NIMBLE discrepancy function is built per
#' specification. The function that comes back then only has to loop over the
#' draws.
#'
#' For each draw it puts the parameter values into the model, puts the dataset
#' into the data nodes, evaluates every discrepancy, simulates a replicate,
#' and evaluates every discrepancy again. The model is left as it was found.
#'
#' The model given here is the one the returned function reads and writes. It
#' does not have to be the model the MCMC runs on, and when running several
#' calibration worlds at once each one should have its own.
#'
#' @param model A NIMBLE model.
#' @param discrepancies One `cppp_discrepancy` object or a list of them.
#' @param simulation A `cppp_simulation` object.
#' @param paramNodes Character vector of the model nodes set from each draw.
#'
#' @return A function `function(MCMCSamples, targetData, control, ...)` that
#'   returns a list with `obs` and `sim`. Both have one row per draw and one
#'   column per discrepancy, with the discrepancy names as column names.
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
      model$calculate()
    })

    for (j in seq_len(nDraws)) {
      values(model, paramNodes)        <- MCMCSamples[j, paramCols]
      values(model, simSpec$dataNodes) <- targetData
      model$calculate()

      for (k in seq_len(K)) obs[j, k] <- discFuns[[k]]$run()

      model$simulate(simSpec$simulateNodes, includeData = TRUE)
      model$calculate()

      for (k in seq_len(K)) rep[j, k] <- discFuns[[k]]$run()
    }

    list(obs = obs, sim = rep)
  }
}
