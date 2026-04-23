
calcDiscFunction_nimble <- nimbleFunction(
  setup = function(MCMCSamples, targetData, control) {
    
    # get dimensions and data
    nSim <- nrow(MCMCSamples)
    nDisc <- length(control$discType)
    discTypes <- control$discType
    
    # get node names
    dataNodes     <- control$dataNodes
    paramNodes    <- control$paramNodes
    expectedNodes <- control$expectedNodes
    
    # get node indices
    paramNodesIndices <- match(paramNodes, colnames(MCMCSamples))
    
    # get nodes to simulate
    nodesSim <- control$model$getDependencies(control$paramNames, self = FALSE)
    
    # create a model copy and get original
    model_copy <- control$model$newModel()
    model <- control$model
    
    # put the helper function generators into the setup. 
    meanFun   <- calc_mean_nimble(model_copy)
    devFun    <- calc_deviance_nimble(model_copy, model)
    chisqFun  <- calc_chisq_nimble(model_copy)
    
  },
  
  run = function() {
    
    # use a 3D array: [Samples, DiscrepancyType, 2]
    # 2 columns for "sim" and "obs"
    discStore <- array(0, dim = c(nSim, nDisc, 2))
    
    for (i in 1:nSim) {
      
      # 1. Update values in both original model and model copy
      values(model_copy, paramNodes) <<- MCMCSamples[i, paramNodesIndices]
      values(model, paramNodes) <<- MCMCSamples[i, paramNodesIndices]
      
      # 2. Simulate
      model_copy$simulate(nodesSim, includeData = TRUE)
      model$simulate(nodesSim, includeData = FALSE)
      
      # 3. Calculate Discrepancies
      for (j in 1:nDisc) {
        if (discTypes[j] == "mean") {
          
          discStore[i, j, ] <- meanFun$run(targetData, dataNodes) 
          
        } else if (discTypes[j] == "deviance") {
          
          discStore[i, j, ] <- devFun$run(dataNodes)
          
        } else if (discTypes[j] == "chisq") {
          
          discStore[i, j, ] <- chisqFun$run(targetData, 
                                            dataNodes, expectedNodes)
        }
      }
    }
    
    returnType(double(3)) 
    return(discStore)
  }
)

calc_mean_nimble <- nimbleFunction(
  
  setup = function(model_copy) {
  },
  
  run = function(targetData = double(1), 
                 dataNodes = character(1)) {
    
    # replicated discrepancies
    discRep <- mean(values(model_copy, dataNodes))
    
    # observed discrepancies
    discObs <- mean(targetData)
    
    returnType(double(1))
    return(c(discRep, discObs))
  }
)

calc_deviance_nimble <- nimbleFunction(
  
  setup = function(model_copy, model) {
  },
  
  run = function(dataNodes = character(1)) {
  
  # replicated discrepancies
  discRep <- 2 * model_copy$calculate(dataNodes)
  
  # observed discrepancies
  discObs <- 2 * model$calculate(dataNodes)
  
  returnType(double(1))
  return(c(discRep, discObs))
  }
)


calc_chisq_nimble <- nimbleFunction(
  
  setup = function(model_copy) {
  },
  
  run = function(targetData = double(1), 
                 dataNodes = character(1), expectedNodes = character(1)) {
  
  # get expected data
  ## AK: check that this works for multiple expected nodes
  data_exp <- values(model_copy, expectedNodes)
  
  data_sim <- values(model_copy, dataNodes)
  
  # replicated discrepancies
  discRep <- sum(
    (data_sim - data_exp) ^ 2 /
      (data_exp + 1e-6)
  )
  
  # observed discrepancies
  discObs <- sum(
    (targetData - data_exp) ^ 2 /
      (data_exp + 1e-6)
  )
  
  returnType(double(1))
  return(c(discRep, discObs))
  }
)

