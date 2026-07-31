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

  K <- length(discs)

  ## The loop over draws lives inside a nimbleFunction, so it runs in compiled
  ## code rather than crossing back into R for every draw and every
  ## discrepancy. Compiling it also compiles the discrepancies with it, in one
  ## call.
  calcNF <- discrepancyCalculatorNF(
    model      = model,
    discs      = discs,
    dataNodes  = simSpec$dataNodes,
    simNodes   = simSpec$simulateNodes,
    paramNodes = paramNodes
  )

  if (compile) {
    compiled <- compileNimble(list(model, calcNF))
    calcNF   <- compiled[[2]]
  }

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

    if (!is.numeric(targetData) || length(targetData) != length(simSpec$dataNodes)) {
      stop("`targetData` must be a numeric vector with one value per data node (",
           length(simSpec$dataNodes), " expected).", call. = FALSE)
    }

    res <- calcNF$run(draws, as.numeric(targetData))

    obs <- matrix(res[, , 1], ncol = K)
    sim <- matrix(res[, , 2], ncol = K)
    colnames(obs) <- discNames
    colnames(sim) <- discNames

    list(obs = obs, sim = sim)
  }
}


#' The draws loop, as a nimbleFunction
#'
#' Holds the discrepancies in a `nimbleFunctionList` so they compile together
#' with the loop that uses them. Not called directly; see
#' [makeDiscrepancyCalculator()].
#'
#' @param model A NIMBLE model.
#' @param discs List of completed `cppp_discrepancy` objects.
#' @param dataNodes Nodes the dataset is written into.
#' @param simNodes Nodes resimulated for a replicate.
#' @param paramNodes Nodes set from each draw.
#' @keywords internal
discrepancyCalculatorNF <- nimbleFunction(
  setup = function(model, discs, dataNodes, simNodes, paramNodes) {

    ## SP: we get dependencies of parameters because we call calculate after changing values
    ## of the model parameters. self=FALSE to avoid overriding of deterministic nodes
    ## (e.g. user monitoring sigma (deterministic) instead of log_sigma (stochastic))
    paramDeps <- model$getDependencies(paramNodes, self = FALSE)
    ## After simulating a replicate: the data's own densities and anything below.
    dataDeps  <- model$getDependencies(dataNodes, self = TRUE)
    restoreDeps <- model$topologicallySortNodes(unique(c(paramDeps, dataDeps)))

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

    ## SP: we need to leave the model as we found it. Saving once here is
    ## enough, because every draw overwrites these before reading them.
    savedData   <- values(model, dataNodes)
    savedParams <- values(model, paramNodes)

    for (i in 1:nDraws) {
      values(model, paramNodes) <<- MCMCOutput[i, ]
      values(model, dataNodes)  <<- targetData
      model$calculate(paramDeps)

      for (k in 1:K) results[i, k, 1] <- discList[[k]]$run()

      model$simulate(simNodes, includeData = TRUE)
      model$calculate(dataDeps)

      for (k in 1:K) results[i, k, 2] <- discList[[k]]$run()
    }

    values(model, dataNodes)  <<- savedData
    values(model, paramNodes) <<- savedParams
    model$calculate(restoreDeps)

    returnType(double(3))
    return(results)
  }
)
