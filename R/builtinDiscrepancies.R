############################################################
## Built-in discrepancies
############################################################
##
## Every discrepancy (built-in or user supplied) is a nimbleFunction that
## declares `contains = discrepancyBase`, whose setup takes three arguments,
## in this order:
##
##     setup = function(model, dataNodes, modelNodes)
##
## and whose run() takes no arguments, reads the model's current state, and
## returns one number (`returnType(double(0))`). A discrepancy that does not
## need `modelNodes` simply ignores it; it is passed as character(0).
##
## Because the structure is fixed, `makeDiscrepancyNimbleFun()` can build a built-in
## and a user-supplied discrepancy with the same single call, and everything
## downstream treats them identically.
##
## A discrepancy that needs particular nodes checks for them in its own setup
## code (which is ordinary R), so the requirement lives with the discrepancy
## instead of in a separate table.
##
## `discrepancyBase` is defined here, above the built-ins that use it, so that
## nothing depends on the order in which R loads the package's files.


#' The shape every discrepancy has
#'
#' A discrepancy is a nimbleFunction whose `run()` takes no arguments and
#' returns a single number. Declaring `contains = discrepancyBase` is what lets
#' the package hold several discrepancies together and compile them in one go.
#'
#' Write your own like this:
#'
#' ```
#' myDisc <- nimbleFunction(
#'   contains = discrepancyBase,
#'   setup = function(model, dataNodes, modelNodes) {},
#'   run = function() {
#'     returnType(double(0))
#'     return(min(values(model, dataNodes)))
#'   }
#' )
#' ```
#'
#' Then name it with [discrepancy()]: `discrepancy("min", fun = myDisc)`.
#'
#' @return A nimbleFunctionVirtual, used as the `contains` argument.
#' @seealso [discrepancy()], [discrepancyBuiltins]
#' @export
discrepancyBase <- nimbleFunctionVirtual(
  run = function() {
    returnType(double(0))
  }
)

#' @keywords internal
meanDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {},
  run = function() {
    returnType(double(0))
    return(mean(values(model, dataNodes)))
  }
)

#' @keywords internal
varianceDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {},
  run = function() {
    returnType(double(0))
    return(var(values(model, dataNodes)))
  }
)

#' @keywords internal
devianceDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {},
  run = function() {
    returnType(double(0))
    ## Deviance = -2 * log-likelihood; calculate() returns the log-likelihood.
    return(-2 * model$calculate(dataNodes))
  }
)

#' @keywords internal
chisquaredDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {
    if (length(modelNodes) != length(dataNodes)) {
      stop("Discrepancy 'chisquared' requires `modelNodes` matching `dataNodes` in length.",
           call. = FALSE)
    }
  },
  run = function() {
    dataVal  <- values(model, dataNodes)
    modelVal <- values(model, modelNodes)
    returnType(double(0))
    return(sum((dataVal - modelVal)^2 / (modelVal + 1e-6)))
  }
)

#' @keywords internal
freemantukeyDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {
    if (length(modelNodes) != length(dataNodes)) {
      stop("Discrepancy 'freemantukey' requires `modelNodes` matching `dataNodes` in length.",
           call. = FALSE)
    }
  },
  run = function() {
    dataVal  <- values(model, dataNodes)
    modelVal <- values(model, modelNodes)
    returnType(double(0))
    return(sum((sqrt(dataVal) - sqrt(modelVal))^2))
  }
)


#' Built-in discrepancies, by name
#'
#' The discrepancies the package ships with. `discrepancy(name)` looks a name up
#' here; `discrepancy(name, fun = )` supplies one of the same shape instead.
#' The list lives in the same file as the functions so that nothing depends on
#' the order in which R loads the package's files.
#'
#' @keywords internal
discrepancyBuiltins <- list(
  mean         = meanDisc,
  variance     = varianceDisc,
  deviance     = devianceDisc,
  chisquared   = chisquaredDisc,
  freemantukey = freemantukeyDisc
)
