#' Create a simulation specification
#'
#' A simulation specification describes how replicated data should be generated
#' from a model when evaluating discrepancies.
#'
#' @param mode Character scalar naming the simulation mode. Current options are
#'   `"conditional"` and `"marginal"`.
#' @param dataNodes Optional character vector of data node names. If `NULL`,
#'   these may be inferred from a NIMBLE model later.
#' @param simulateNodes Optional character vector of nodes to resimulate when
#'   generating a replicate. If `NULL`, defaults depend on `mode`.
#'
#' @return An object of class `cppp_simulation`.
#' @export
simulation <- function(mode = "conditional",
                       dataNodes = NULL,
                       simulateNodes = NULL) {
  validModes <- c("conditional", "marginal")

  if (!is.character(mode) || length(mode) != 1L) {
    stop("`mode` must be a single character string.", call. = FALSE)
  }

  if (!(mode %in% validModes)) {
    stop(
      sprintf(
        "`mode` must be one of: %s.",
        paste(validModes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  x <- list(
    mode = mode,
    dataNodes = dataNodes,
    simulateNodes = simulateNodes
  )

  class(x) <- "cppp_simulation"
  x
}


#' Complete a simulation specification from a model
#'
#' Fill in defaults for a simulation specification using information from a
#' NIMBLE model.
#'
#' @param model A NIMBLE model.
#' @param sim A `cppp_simulation` object.
#'
#' @return A completed `cppp_simulation` object.
#' @keywords internal
completeSimulation <- function(model, sim) {
  stopifnot(inherits(sim, "cppp_simulation"))

  out <- sim

  if (is.null(out$dataNodes)) {
    out$dataNodes <- model$getNodeNames(dataOnly = TRUE)
  }

  out$dataNodes <- model$expandNodeNames(
    out$dataNodes,
    returnScalarComponents = TRUE
  )

  if (is.null(out$simulateNodes)) {
    if (out$mode == "conditional") {
      out$simulateNodes <- out$dataNodes
    } else if (out$mode == "marginal") {
      stop(
          "For `mode = \"marginal\"`, `simulateNodes` must be supplied.",
        call. = FALSE
      )
    }
  }

  out$simulateNodes <- model$expandNodeNames(
    out$simulateNodes,
    returnScalarComponents = TRUE
  )

  out
}
