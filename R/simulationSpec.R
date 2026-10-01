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
#' @param latentNodes Character vector of the latent nodes to redraw, needed
#'   for `"marginal"` and not used for `"conditional"`. Name the top of what to
#'   redraw: everything below these nodes, data included, is redrawn too.
#'
#' @return An object of class `cppp_simulation`.
#' @export
simulation <- function(mode = "conditional",
                       latentNodes = NULL) {
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

  if (mode == "conditional" && !is.null(latentNodes)) {
    stop("`latentNodes` is only used with `mode = \"marginal\"`.", call. = FALSE)
  }

  x <- list(
    mode = mode,
    latentNodes = latentNodes
  )

  class(x) <- "cppp_simulation"
  x
}


#' Work out the node lists from a model
#'
#' This function groups the model nodes the packagage needs for the entire cppp
#' procedure. The calculator, the simulate function and the MCMC all take
#' their nodes from here.
#'
#' @param model A NIMBLE model.
#' @param dataNodes Optional character vector of the observed nodes. If `NULL`,
#'   all data nodes in the model.
#' @param paramNodes Optional character vector of the nodes set from each
#'   posterior draw. If `NULL`, all stochastic nodes that are not data.
#' @param simulation A [simulation()] specification. If `NULL`,
#'   `simulation("conditional")`.
#'
#' @return A list of class `cppp_nodes` with
#' * `data`: the  nodes, one entry per value.
#' * `params`: the nodes set from each draw, one entry per value.
#' * `simulate`: the nodes redrawn for a replicate, in model order.
#' * `paramDeps`: the nodes to recalculate after setting a draw.
#' * `saved`: every node a replicate can change, which the calculator puts
#'   back when it is done.
#' @keywords internal
completeNodes <- function(model, dataNodes = NULL, paramNodes = NULL,
                          simulation = NULL) {

  if (is.null(simulation)) simulation <- simulation("conditional")
  stopifnot(inherits(simulation, "cppp_simulation"))

  if (is.null(dataNodes)) dataNodes <- model$getNodeNames(dataOnly = TRUE)
  data <- model$expandNodeNames(dataNodes, returnScalarComponents = TRUE)
  if (!all(data %in% model$getNodeNames(stochOnly = TRUE))) {
    stop("All data nodes must be stochastic nodes in the model.", call. = FALSE)
  }

  if (is.null(paramNodes)) {
    paramNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  }
  params <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
  if (length(params) == 0L) {
    stop("`paramNodes` did not match any node in the model.", call. = FALSE)
  }

  if (simulation$mode == "conditional") {
    simulate <- data
  } else {
    ## Only the user knows which nodes are the latent states.
    if (is.null(simulation$latentNodes)) {
      stop("For `mode = \"marginal\"`, give `latentNodes`: the latent nodes to redraw.",
           call. = FALSE)
    }
    latent <- model$expandNodeNames(simulation$latentNodes,
                                    returnScalarComponents = TRUE)
    if (any(latent %in% data)) {
      stop("`latentNodes` must not include data nodes: ",
           paste(intersect(latent, data), collapse = ", "), call. = FALSE)
    }
    ## Simulating a node leaves the nodes below it alone, so redraw everything
    ## below the latents too, including the helper nodes NIMBLE adds and the
    ## data. `downstream = TRUE` goes past the first stochastic node it meets.
    simulate <- model$getDependencies(c(latent, data), self = TRUE,
                                      downstream = TRUE)
  }

  ## self = FALSE so a derived parameter (sigma monitored instead of
  ## log_sigma) keeps the value written from the draw.
  paramDeps <- model$getDependencies(params, self = FALSE)

  out <- list(
    data      = data,
    params    = params,
    simulate  = simulate,
    paramDeps = paramDeps,
    saved     = unique(c(params, paramDeps, simulate))
  )
  class(out) <- "cppp_nodes"
  out
}
