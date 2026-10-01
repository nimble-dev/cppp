#' Build a discrepancy calculator
#'
#' Create a function that computes observed and replicated discrepancy values
#' for a set of posterior draws, from discrepancy and simulation
#' specifications.
#'
#' @details
#' For each posterior draw, the returned function sets the parameter nodes,
#' evaluates every discrepancy using the data in the model,
#' simulates a replicate dataset, and evaluates every discrepancy
#' again. The model is restored to its original state on completion.
#'
#' The returned function is suitable as the `discFun` argument of
#' [runCalibration()] and [runCalibrationNIMBLE()]. It may also be called
#' directly, to inspect discrepancy values before running a calibration.
#'
#' `model` must be uncompiled, and should be an instance not in use elsewhere,
#' such as `model$newModel()`, since the calculator writes parameter values and
#' data into it.
#'
#' @param model A NIMBLE model.
#' @param discrepancies A [discrepancy()] specification, or a list of them.
#' @param simulation A [simulation()] specification.
#' @param paramNodes Character vector naming the model nodes to set from each
#'   posterior draw. These must appear among the column names of the draws.
#' @param dataNodes Optional character vector of the data nodes the dataset is
#'   written into. If `NULL`, all data nodes in the model.
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
                                      dataNodes = NULL, compile = TRUE) {

  nodes <- completeNodes(model, dataNodes, paramNodes, simulation)
  parts <- buildDiscrepancyCalculator(model, discrepancies, nodes)

  if (compile) {
    compiled     <- compileNimble(list(model, parts$calcNF))
    parts$calcNF <- compiled[[2]]
  }

  wrapDiscrepancyCalculator(parts)
}


#' Build the calculator without compiling
#'
#' The first half of [makeDiscrepancyCalculator()]. It checks the
#' specifications and builds the nimbleFunction that holds the loop over draws,
#' but does not compile it.
#'
#' Kept separate so the nimbleFunction can be compiled together with the model
#' and the MCMC in one call, rather than in a call of its own. Hand the result
#' to [wrapDiscrepancyCalculator()], with `calcNF` replaced by the compiled
#' version if there is one.
#'
#' @param model An uncompiled NIMBLE model.
#' @param discrepancies A [discrepancy()] specification, or a list of them.
#' @param nodes The node lists, from [completeNodes()].
#'
#' @return A list with the nimbleFunction `calcNF`, the expanded `paramNodes`,
#'   the `dataNodes` the dataset is written into, and the `discNames`.
#' @keywords internal
buildDiscrepancyCalculator <- function(model, discrepancies, nodes) {

  ## A nimbleFunction is built against an uncompiled model, so we need one to
  ## start from even though the work then happens on the compiled copy.
  if (inherits(model, "CmodelBaseClass")) {
    stop("`model` must be an uncompiled NIMBLE model; the calculator compiles its own copy.",
         call. = FALSE)
  }

  discs   <- standardizeDiscrepancies(discrepancies)
  discs   <- lapply(discs, function(d) completeDiscrepancy(model, d))

  discNames <- vapply(discs, function(d) d$name, character(1))
  if (anyDuplicated(discNames)) {
    stop("Each discrepancy needs its own name; repeated: ",
         paste(unique(discNames[duplicated(discNames)]), collapse = ", "),
         call. = FALSE)
  }

  ## SP: a discrepancy may look at part of the data only, but never at nodes
  ## outside the data we write into and simulate.
  for (d in discs) {
    unwritten <- setdiff(d$dataNodes, nodes$data)
    if (length(unwritten) > 0L) {
      stop("Discrepancy '", d$name, "' reads data nodes the simulation does not set: ",
           paste(unwritten, collapse = ", "),
           ". Give `discrepancy()` only nodes among the data nodes.",
           call. = FALSE)
    }
  }

  ## We loop over posterior draws withing the nimbleFunction, so it runs in compiled
  ## code rather than crossing back into R for every draw and every
  ## discrepancy. Compiling it also compiles the discrepancies with it, in one
  ## call.
  calcNF <- discrepancyCalculatorNF(
    model = model,
    discs = discs,
    nodes = nodes
  )

  list(
    calcNF     = calcNF,
    paramNodes = nodes$params,
    dataNodes  = nodes$data,
    discNames  = discNames
  )
}


#' Wrap the calculator's pieces into the function you call
#'
#' The second half of [makeDiscrepancyCalculator()]. Takes what
#' [buildDiscrepancyCalculator()] produced and returns the function that
#' [runCalibration()] calls, which lines the draws up for the nimbleFunction
#' and puts names back on the results.
#'
#' Works whether `calcNF` has been compiled or not.
#'
#' @param parts The list returned by [buildDiscrepancyCalculator()].
#'
#' @return A function of `(MCMCSamples, targetData, control, ...)`, returning a
#'   list of two matrices, `obs` and `sim`.
#' @keywords internal
wrapDiscrepancyCalculator <- function(parts) {

  calcNF     <- parts$calcNF
  paramNodes <- parts$paramNodes
  dataNodes  <- parts$dataNodes
  discNames  <- parts$discNames
  K          <- length(discNames)

  function(MCMCSamples, targetData, control = NULL, ...) {

    MCMCSamples <- as.matrix(MCMCSamples)

    ## Compiled code cannot pick columns out by name, so line them up here and
    ## hand over just the parameter columns, in `paramNodes` order.
    absent <- setdiff(paramNodes, colnames(MCMCSamples))
    if (length(absent) > 0L) {
      stop("`MCMCSamples` has no column for: ", paste(absent, collapse = ", "),
           call. = FALSE)
    }
    draws <- MCMCSamples[, match(paramNodes, colnames(MCMCSamples)), drop = FALSE]
    storage.mode(draws) <- "double"

    if (!is.numeric(targetData) || length(targetData) != length(dataNodes)) {
      stop("`targetData` must be a numeric vector with one value per data node (",
           length(dataNodes), " expected).", call. = FALSE)
    }

    res <- calcNF$run(draws, as.numeric(targetData))

    obs <- matrix(res[, , 1], ncol = K)
    sim <- matrix(res[, , 2], ncol = K)
    colnames(obs) <- discNames
    colnames(sim) <- discNames

    list(obs = obs, sim = sim)
  }
}



#' Compute discrepancies for every posterior draw
#'
#' The work behind [makeDiscrepancyCalculator()]. For each draw it sets the
#' parameters, evaluates every discrepancy on the dataset, simulates a
#' replicate, and evaluates them again.
#'
#' Written as a nimbleFunction so the loop runs in compiled code, and so the
#' discrepancies (held in a `nimbleFunctionList`) compile together with it in
#' a single call. Not called directly.
#'
#' @param model A NIMBLE model.
#' @param discs List of completed `cppp_discrepancy` objects.
#' @param nodes The node lists, from [completeNodes()].
#' @keywords internal
discrepancyCalculatorNF <- nimbleFunction(
  setup = function(model, discs, nodes) {

    paramNodes <- nodes$params
    dataNodes  <- nodes$data
    simNodes   <- nodes$simulate
    paramDeps  <- nodes$paramDeps

    ## mvSaved contains the state of the model (nodes values and logProbs etc. )
    ## Only the nodes changes in the run part are saved. These are parameters,
    ## their dependencies, and what we resimulate (see completeNodes()).
    ## Restoring is a copy back rather than a recalculation.
    savedNodes <- nodes$saved
    mvSaved    <- modelValues(model)

    discList <- nimbleFunctionList(discrepancyBase)
    for (i in seq_along(discs)) {
      discList[[i]] <- makeDiscrepancyNimbleFun(model, discs[[i]])
    }
    K <- length(discs)
  },
  run = function(MCMCOutput = double(2), targetData = double(1)) {
    nDraws <- dim(MCMCOutput)[1]

    ## results[draw, discrepancy, 1] is the observed side, [, , 2] the
    ## replicated one.
    results <- array(0, c(nDraws, K, 2))

    ## SP: we need to leave the model as we found it before running the calculations
    nimCopy(from = model, to = mvSaved, row = 1, nodes = savedNodes, logProb = TRUE)

    for (i in 1:nDraws) {
      values(model, paramNodes) <<- MCMCOutput[i, ]
      values(model, dataNodes)  <<- targetData
      model$calculate(paramDeps)

      for (k in 1:K) results[i, k, 1] <- discList[[k]]$run()

      model$simulate(simNodes, includeData = TRUE)
      model$calculate(simNodes)

      for (k in 1:K) results[i, k, 2] <- discList[[k]]$run()
    }

    nimCopy(from = mvSaved, to = model, row = 1, nodes = savedNodes, logProb = TRUE)

    returnType(double(3))
    return(results)
  }
)
