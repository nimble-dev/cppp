#' Create a discrepancy specification
#'
#' A discrepancy specification is simply a list containing the description of a
#' discrepancy. The package can use it to complete missing defaults from a model
#' and build a NIMBLE implementation from the completed specification.
#'
#' There is one way to define a discrepancy: a nimbleFunction whose setup takes
#' `(model, dataNodes, modelNodes)` and whose `run()` reads the model's current
#' state and returns a single number. The package ships a few such functions
#' (see [discrepancyBuiltins]); naming one of them is a shortcut for supplying
#' it yourself. Pass `fun` to use your own.
#'
#' To use an R function inside your own discrepancy, wrap it with
#' `nimble::nimbleRcall()` in your own script, for example:
#'
#' ```
#' sortR <- nimbleRcall(prototype  = function(x = double(1)){},
#'                      Rfun       = "sort",
#'                      returnType = double(1))
#' ```
#'
#' @param name Character scalar naming the discrepancy. It labels this
#'   discrepancy in the results, so it must be supplied even when `fun` is.
#'   Without `fun`, it must name one of the built-in discrepancies.
#' @param dataNodes Optional character vector of data node names. If `NULL`,
#'   these may be inferred from a NIMBLE model.
#' @param modelNodes Optional character vector of model-based nodes used by the
#'   discrepancy besides the data nodes. For example, expected-value nodes for
#'   chi-squared.
#' @param fun Optional nimbleFunction implementing your own discrepancy, in the
#'   shape described above. If `NULL`, `name` is looked up among the built-ins.
#'
#' @return An object of class `cppp_discrepancy`.
#' @export
discrepancy <- function(name,
                        dataNodes = NULL,
                        modelNodes = NULL,
                        fun = NULL) {

  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single character string.", call. = FALSE)
  }

  if (is.null(fun)) {
    if (!(name %in% names(discrepancyBuiltins))) {
      stop(
        sprintf(
          "Unknown discrepancy '%s'. Built-in discrepancies are: %s. Supply `fun` to use your own.",
          name,
          paste(names(discrepancyBuiltins), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  } else {
    if (!is.function(fun)) {
      stop("`fun` must be a nimbleFunction.", call. = FALSE)
    }
    ## Otherwise `name` would label the results with a built-in's name while
    ## computing something else.
    if (name %in% names(discrepancyBuiltins)) {
      stop(
        sprintf(
          "'%s' is a built-in discrepancy; choose a different `name` for your own.",
          name
        ),
        call. = FALSE
      )
    }
  }

  x <- list(
    name = name,
    dataNodes = dataNodes,
    modelNodes = modelNodes,
    fun = fun
  )

  class(x) <- "cppp_discrepancy"
  x
}


#' Standardize discrepancy input
#'
#' Accept either one discrepancy specification or a list of them and always
#' return a list of discrepancy specifications.
#'
#' @param x A `cppp_discrepancy` object or a list of such objects.
#'
#' @return A list of `cppp_discrepancy` objects.
#' @keywords internal
standardizeDiscrepancies <- function(x) {
  if (inherits(x, "cppp_discrepancy")) {
    return(list(x))
  }

  if (is.list(x) &&
      length(x) > 0L &&
      all(vapply(x, inherits, logical(1), "cppp_discrepancy"))) {
    return(x)
  }

  stop(
    "`discrepancy` must be a discrepancy object or a list of discrepancy objects.",
    call. = FALSE
  )
}


#' Complete a discrepancy specification from a model
#'
#' Fill in the node names a discrepancy needs, using information from a NIMBLE
#' model. Nothing here depends on *which* discrepancy it is: a discrepancy that
#' needs particular nodes checks for them itself, when it is built.
#'
#' @param model A NIMBLE model.
#' @param disc A `cppp_discrepancy` object.
#'
#' @return A completed `cppp_discrepancy` object.
#' @keywords internal
completeDiscrepancy <- function(model, disc) {
  stopifnot(inherits(disc, "cppp_discrepancy"))

  out <- disc

  if (is.null(out$dataNodes)) {
    out$dataNodes <- model$getNodeNames(dataOnly = TRUE)
  }

  out$dataNodes <- model$expandNodeNames(
    out$dataNodes,
    returnScalarComponents = TRUE
  )

  ## character(0) rather than NULL: every discrepancy's setup is called with all
  ## three arguments, and those that ignore `modelNodes` must still receive it.
  out$modelNodes <- if (is.null(out$modelNodes)) {
    character(0)
  } else {
    model$expandNodeNames(out$modelNodes, returnScalarComponents = TRUE)
  }

  out
}


#' Build a NIMBLE discrepancy evaluator
#'
#' Create the NIMBLE function for a discrepancy specification, pointed at the
#' given model. The result has a `run()` method that evaluates the discrepancy
#' for the model's current state.
#'
#' Built-in and user-supplied discrepancies are built by the same call; the only
#' difference is where the nimbleFunction came from.
#'
#' @param model A NIMBLE model.
#' @param disc A `cppp_discrepancy` object.
#'
#' @return A NIMBLE function object with a scalar `run()` method.
#' @export
makeDiscrepancyNimbleFun <- function(model, disc) {
  disc <- completeDiscrepancy(model, disc)

  gen <- if (!is.null(disc$fun)) {
    disc$fun
  } else {
    discrepancyBuiltins[[disc$name]]
  }

  if (is.null(gen)) {
    stop(
      sprintf(
        "Unknown discrepancy '%s'. Built-in discrepancies are: %s.",
        disc$name,
        paste(names(discrepancyBuiltins), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gen(model, disc$dataNodes, disc$modelNodes)
}
