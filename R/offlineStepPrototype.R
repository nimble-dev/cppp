#' Run one offline discrepancy step in R
#'
#' This is a small prototype for one iteration of an offline calibration loop.
#' It assumes the model is already in the state of one posterior draw, or
#' optionally allows selected model nodes to be updated from `thetaRow` first.
#'
#' @param model A NIMBLE model.
#' @param disc A `cppp_discrepancy` object.
#' @param sim A `cppp_simulation` object.
#' @param thetaRow Optional numeric vector giving one posterior draw.
#' @param paramNodes Optional character vector of model nodes updated from
#'   `thetaRow`. Required if `thetaRow` is supplied.
#'
#' @return A list with components `obs` and `rep`.
#' @export
runOneDiscrepancyStepR <- function(model,
                                   disc,
                                   sim,
                                   thetaRow = NULL,
                                   paramNodes = NULL) {
  disc <- completeDiscrepancy(model, disc)
  sim <- completeSimulation(model, sim)
  ## The NIMBLE discrepancy is the single source of truth; run it uncompiled
  ## here. It reads the model's current state, so we set the state first and
  ## then call its run() method.
  discFun <- makeDiscrepancyNimbleFun(model, disc)

  savedSimState <- values(model, sim$simulateNodes)
  savedParamState <- NULL

  on.exit({
    values(model, sim$simulateNodes) <- savedSimState
    if (!is.null(paramNodes) && !is.null(savedParamState)) {
      values(model, paramNodes) <- savedParamState
    }
    model$calculate()
  }, add = TRUE)

  if (!is.null(thetaRow)) {
    if (is.null(paramNodes)) {
      stop("`paramNodes` must be supplied if `thetaRow` is supplied.", call. = FALSE)
    }

    paramNodes <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    savedParamState <- values(model, paramNodes)

    if (is.null(names(thetaRow))) {
      if (length(thetaRow) != length(paramNodes)) {
        stop("Unnamed `thetaRow` must have the same length as `paramNodes`.", call. = FALSE)
      }
      values(model, paramNodes) <- thetaRow
    } else {
      if (!all(paramNodes %in% names(thetaRow))) {
        stop("Named `thetaRow` must contain all `paramNodes`.", call. = FALSE)
      }
      values(model, paramNodes) <- thetaRow[paramNodes]
    }
    model$calculate()
  }

  obs <- discFun$run()

  model$simulate(sim$simulateNodes, includeData = TRUE)
  model$calculate()
  rep <- discFun$run()

  list(obs = obs, rep = rep)
}
