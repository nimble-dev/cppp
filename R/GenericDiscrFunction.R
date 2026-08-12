calcDiscrepancies <- nimbleFunction(
  setup = function(control) {

    # get model
    model <- control$model

    # get dimensions and data
    nDisc <- length(control$discType)
    discTypes <- control$discType

    # get node names
    dataNodes     <- control$dataNodes
    paramNodes    <- control$paramNodes
    expectedNodes <- control$expectedNodes

    # get parameter and data dependencies
    paramDependencies <- model$getDependencies(paramNodes)
    dataDependencies  <- model$getDependencies(dataNodes)

    # get node indices
    paramNodeIndices <- match(paramNodes, colnames(MCMCSamples))

    # get nodes to simulate
    simNodes <- model$getDependencies(paramNodes, self = FALSE)
    simNodes <- model$topologicallySortNodes(simNodes)

    # create empty disc function list
    discrepancyFunctionList <- nimbleFunctionList(discrepancyFunction_BASE)

    # use look-up table to add disc functions
    for (i in seq_along(discTypes)) {

      if (discTypes[i] %in% names(discLookup)) {
        # grab the function from the lookup and initialize it
        discrepancyFunctionList[[i]] <- discLookup[[discTypes[i]]](model,
                                                                   control)
      } else {
        stop(paste0("discType '", discTypes[i],
                    "' is invalid. Must be: mean, variance, deviance, ",
                    "chisquared, or freemantukey"))
      }

    }

  },
  run = function(MCMCSamples = double(2)) {

    origDataValues  <- values(model, dataNodes)
    origParamValues <- values(model, paramNodes)
    nSamples <- dim(MCMCSamples)[1]

    ## results[j, i, k]:
    ## j = discrepancy index,
    ## i = MCMC iteration,
    ## k = 1 for observed or 2 for simulated
    results <- array(dim = c(nDisc, nSamples, 2) )

    for (i in 1:nSamples) {
      ## put MCMC values from sampled iteration in the model
      values(model, paramNodes) <<- MCMCSamples[i, paramNodeIndices]

      ## calculate
      model$calculate(paramDependencies)

      ## calculate observed discrepancies
      for (j in 1:nDisc) {
        obsDiscrepancy <- discrepancyFunctionList[[j]]$run()
        results[j, i, 1] <- obsDiscrepancy
      }

      ## simulate from posterior predictive
      model$simulate(simNodes, includeData = TRUE)
      model$calculate(dataDependencies)

      for (j in 1:nDisc) {
        simDiscrepancy <- discrepancyFunctionList[[j]]$run()
        results[j, i, 2] <- simDiscrepancy
      }

      ## Careful here with latent states!
      ## Return the model to original state (i.e. state upon entry to this function)
      values(model, dataNodes) <<- origDataValues
      model$calculate(dataDependencies)
    }
    values(model, paramNodes) <<- origParamValues
    model$calculate(paramDependencies)

    returnType(double(3))
    return(results)
  }
)
