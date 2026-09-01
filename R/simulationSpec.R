#' Create a simulation specification
#'
#' A simulation specification is a list that describes how replicated data should
#' be generated from a model when calculating discrepancies.
#'
#' The are two modes:
#' * `"conditional"` latent quantities are fixed at the values from the
#' posterior draw, so that only data is resimulated
#'(indeed conditional to latent qunatities).
#' * `"marginal"` also redraws the latent quantities, so the replicate
#' comes from the model integrated over them.
#'
#' @param mode Character scalar specifying the the mode (`"conditional"` or `"marginal"`)
#' @param dataNodes Optional character vector of data node names. If `NULL`,
#'   these may be inferred from a NIMBLE model.
#' @param simulateNodes Optional character vector of nodes to resimulate when
#'   generating a replicate. If `NULL`, defaults depend on `mode`: the data
#'   nodes for `"conditional"`, and every stochastic node downstream of the
#'   parameters for `"marginal"`.
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
      ## Everything that is not data is set equal to the posterior draw.
      out$simulateNodes <- out$dataNodes
    } else if (out$mode == "marginal") {
      ## SP: no default. User must provide latent states
      stop(
        "For `mode = \"marginal\"`, `simulateNodes` must be supplied: list the latent nodes to redraw as well as the data nodes.",
        call. = FALSE
      )
    }
  }

  out$simulateNodes <- model$expandNodeNames(
    out$simulateNodes,
    returnScalarComponents = TRUE
  )

  out$simulateNodes <- model$topologicallySortNodes(out$simulateNodes)

  out
}
