calcDiscFunction <- function(MCMCSamples  = MCMCSamples,
                             targetData   = observedData,
                             control      = discControl) {
  
  # create empty list of matrices
  disc <- setNames(replicate(length(discType), 
                             matrix(NA, nrow(MCMCSamples), 2, 
                                    dimnames = list(NULL, 
                                                    c("sim", "obs"))), 
                             simplify = FALSE), 
                   discType)
  
  # make a copy of model for simulating replicated discrepancies
  model_copy <- control$model$newModel()
  
  # get nodes to simulate
  nodesSim <- model_copy$getDependencies(control$paramNames, self = FALSE)
  
  # loop through MCMCsamples
  for (i in 1:nrow(MCMCSamples)) {
    
    ##
    # simulate nodes for replicated discrepancies
    ##
    
    # add values to model copy
    values(model_copy, control$paramNodes) <- MCMCSamples[i, control$paramNodes]
    
    # simulate other nodes
    model_copy$simulate(nodesSim, includeData = TRUE)
    
    ##
    # simulate nodes for observed discrepancies
    ##
    
    # add values to model
    values(control$model, control$paramNodes) <- MCMCSamples[i, control$paramNodes]
    
    # simulate other nodes (not including data)
    control$model$simulate(nodesSim, includeData = FALSE)
    
    # loop through discrepancy types
    for (j in seq_along(control$discType)) {
      
      disc[[j]][i, ] <- switch(
        control$discType[j],
        "mean"     = calc_mean(model_copy, targetData, control),
        "deviance" = calc_deviance(model_copy, control),
        "chisq"    = calc_chisq(model_copy, targetData, control),
        # AK: should probably put this earlier
        # Default case (the error)
        stop(paste0("discType ", control$discType[j], 
                    " is invalid. Must be: mean, deviance, chisq"))
      )
      
    }
    
  }
  
  # AK: does this return type make sense? 
  # List of length(number of discrepancies)
  # Each element of the list is a matrix with nrow = # MCMCSamples, ncol = 2
  return(disc)
}

calc_mean <- function(model_copy, targetData, control) {
  
  # replicated discrepancies
  discRep <- mean(values(model_copy, control$dataNodes))
  
  # observed discrepancies
  discObs <- mean(targetData)
  
  return(c(discRep, discObs))
}

calc_deviance <- function(model_copy, control) {
  
  # replicated discrepancies
  discRep <- 2 * model_copy$calculate(control$dataNodes)
  
  # observed discrepancies
  discObs <- 2 * control$model$calculate(control$dataNodes)
  
  return(c(discRep, discObs))
  
}

calc_chisq <- function(model_copy, targetData, control) {
  
  # get expected data
  ## AK: check that this works for multiple expected nodes
  data_exp <- values(model_copy, control$expectedNodes)
  
  # replicated discrepancies
  discRep <- sum(
    (values(model_copy, control$dataNodes) - data_exp) ^ 2 /
      (data_exp + 1e-6)
  )
  
  # observed discrepancies
  discObs <- sum(
    (targetData - data_exp) ^ 2 /
      (data_exp + 1e-6)
  )
  
  return(c(discRep, discObs))
  
}
