##-------------------------------------------------------##
## This file may contain some already implemented discrepancies
##-------------------------------------------------------##

## discrepancy base class
discrepancyFunction_BASE <- nimbleFunctionVirtual(
  run = function() returnType(double())
)

# mean discrepancy function
MeanDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, control) {

    data <- control$dataNodes

  },
  run = function() {

    disc <- mean(values(model, data))

    returnType(double(0))
    return(disc)
  }
)

# variance discrepancy function
VarDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, control) {

    data <- control$dataNodes

  },
  run = function() {

    disc <- var(values(model, data))

    returnType(double(0))
    return(disc)
  }
)

# deviance discrepancy function
DevianceDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, control) {

    data <- control$dataNodes

  },
  run = function() {

    disc <- 2 * model$calculate(data)

    returnType(double(0))
    return(disc)
  }
)

# chi-squared discrepancy function
ChiSqDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, control) {

    dataNodes <- control$dataNodes
    expectedNodes <- control$expectedNodes

  },
  run = function() {

    data_val <- values(model, dataNodes)

    ## AK: check that this works for multiple expected nodes
    data_exp <- values(model, expectedNodes)

    disc <- sum(
      (data_val - data_exp) ^ 2 /
        (data_exp + 1e-6)
    )

    returnType(double(0))
    return(disc)
  }
)

# Freeman-Tukey discrepancy function
FTukeyFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, control) {

    dataNodes <- control$dataNodes
    expectedNodes <- control$expectedNodes

  },
  run = function() {

    data_val <- values(model, dataNodes)

    ## AK: check that this works for multiple expected nodes
    data_exp <- values(model, expectedNodes)

    disc <- sum(
      (sqrt(data_val) - sqrt(data_exp)) ^ 2
    )

    returnType(double(0))
    return(disc)
  }
)

# discrepancy look-up map
discLookup <- list(
  mean = MeanDiscFunction,
  variance = VarDiscFunction,
  deviance = DevianceDiscFunction,
  chisquared = ChiSqDiscFunction,
  freemantukey = FTukeyFunction
)
