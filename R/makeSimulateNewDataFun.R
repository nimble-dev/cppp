#' Build a function that simulates one replicate dataset
#'
#' Builds the function that makes one replicate dataset from one posterior draw,
#' using the node lists from [completeNodes()].
#'
#' Give it its own copy of the model, `model$newModel()`. It writes parameter
#' values into the model and simulates into the data nodes.
#'
#' @param model A NIMBLE model.
#' @param nodes The node lists, from [completeNodes()].
#'
#' @return A function `function(thetaRow, control = NULL, ...)` returning the
#'   replicate dataset as a numeric vector, one value per data node.
#' @seealso [makeDiscrepancyCalculator()], [runCalibrationNIMBLE()]
#' @keywords internal
makeSimulateNewDataFun <- function(model, nodes) {

  paramNodes <- nodes$params
  paramDeps  <- nodes$paramDeps

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

    model$simulate(nodes$simulate, includeData = TRUE)

    values(model, nodes$data)
  }
}
