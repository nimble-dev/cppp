#' Build a function that simulates one replicate dataset
#'
#' A calibration needs a replicate dataset for each of its worlds: take one
#' posterior draw, put it into the model, and simulate. That is what the
#' [simulation()] specification already describes, so it can be built rather
#' than written by hand.
#'
#' Give it its own copy of the model, `model$newModel()`. It writes parameter
#' values into the model and simulates into the data nodes.
#'
#' @param model A NIMBLE model.
#' @param simulation A [simulation()] specification.
#' @param paramNodes Character vector naming the model nodes set from the draw.
#'
#' @return A function `function(thetaRow, control = NULL, ...)` returning the
#'   replicate dataset as a numeric vector, one value per data node.
#' @seealso [makeDiscrepancyCalculator()], [runCalibrationNIMBLE()]
#' @export
makeSimulateNewDataFun <- function(model, simulation, paramNodes) {

  simSpec    <- completeSimulation(model, simulation)
  paramNodes <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)

  ## Below the parameters, never the parameters themselves: a parameter can be
  ## a derived node, such as sigma when the prior is on log(sigma), and
  ## recalculating one of those would rewrite the value the draw just gave us.
  ## Same reasoning as in makeDiscrepancyCalculator().
  paramDeps <- model$getDependencies(paramNodes, self = FALSE)

  function(thetaRow, control = NULL, ...) {

    if (is.null(names(thetaRow))) {
      if (length(thetaRow) != length(paramNodes)) {
        stop("Unnamed `thetaRow` must have one value per parameter node (",
             length(paramNodes), " expected).", call. = FALSE)
      }
      draw <- thetaRow
    } else {
      absent <- setdiff(paramNodes, names(thetaRow))
      if (length(absent) > 0L) {
        stop("`thetaRow` has no value for: ", paste(absent, collapse = ", "),
             call. = FALSE)
      }
      draw <- thetaRow[paramNodes]
    }

    values(model, paramNodes) <- as.numeric(draw)
    model$calculate(paramDeps)

    model$simulate(simSpec$simulateNodes, includeData = TRUE)

    values(model, simSpec$dataNodes)
  }
}
