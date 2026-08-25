#' Build a function that simulates one replicate dataset
#'
#' Builds the function that makes one replicate dataset from one posterior draw,
#' following a [simulation()] specification.
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
#' @keywords internal
makeSimulateNewDataFun <- function(model, simulation, paramNodes) {

  simSpec    <- completeSimulation(model, simulation)
  paramNodes <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)

  ## Get dependencies of the parameters excluding the parameters to avoid
  ## overwriting of lifted nodes
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
