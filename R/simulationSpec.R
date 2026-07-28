#' Create a simulation specification
#'
#' A simulation specification describes how replicated data should be generated
#' from a model when evaluating discrepancies.
#'
#' The two modes differ in how much of the model is redrawn for a replicate.
#' `"conditional"` holds any latent quantities fixed at the values from the
#' posterior draw and resimulates only the data. `"marginal"` also redraws the
#' latent quantities, so the replicate comes from the model integrated over
#' them.
#'
#' @param mode Character scalar naming the simulation mode. Current options are
#'   `"conditional"` and `"marginal"`.
#' @param dataNodes Optional character vector of data node names. If `NULL`,
#'   these may be inferred from a NIMBLE model later.
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
#' @param paramNodes Optional character vector of the model nodes set from each
#'   posterior draw. Used only in `"marginal"` mode, to work out which nodes to
#'   resimulate. Not needed if `simulateNodes` is already supplied.
#'
#' @return A completed `cppp_simulation` object.
#' @keywords internal
completeSimulation <- function(model, sim, paramNodes = NULL) {
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
      ## Everything above the data stays at the posterior draw.
      out$simulateNodes <- out$dataNodes
    } else if (out$mode == "marginal") {
      if (is.null(paramNodes)) {
        stop(
          "For `mode = \"marginal\"`, supply `paramNodes` so the nodes to resimulate can be worked out, or supply `simulateNodes` yourself.",
          call. = FALSE
        )
      }
      ## Redraw everything the parameters feed into, latent nodes included, not
      ## just the data. getDependencies() returns them in an order it is safe to
      ## simulate in (parents before children).
      out$simulateNodes <- model$getDependencies(
        paramNodes,
        self       = FALSE,
        stochOnly  = TRUE,
        downstream = TRUE,
        includeData = TRUE
      )

      if (length(out$simulateNodes) == 0L) {
        stop(
          "No stochastic nodes depend on `paramNodes`, so there is nothing to resimulate in marginal mode.",
          call. = FALSE
        )
      }
    }
  }

  out$simulateNodes <- model$expandNodeNames(
    out$simulateNodes,
    returnScalarComponents = TRUE
  )

  out
}
