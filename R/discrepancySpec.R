#' Create a discrepancy specification
#'
#' A discrepancy specification is a lightweight user-facing description of a
#' discrepancy. The package can later complete missing defaults from a model
#' and build an R or NIMBLE implementation from the completed specification.
#'
#' @param kind Character scalar naming the structural kind of discrepancy.
#'   Current options are `"data_only"`, `"data_plus_model"`, and
#'   `"data_plus_external"`.
#' @param name Character scalar naming the discrepancy.
#' @param dataNodes Optional character vector of data node names. If `NULL`,
#'   these may be inferred from a NIMBLE model later.
#' @param modelNodes Optional character vector of model-based nodes used by the
#'   discrepancy besides the data nodes. For example, expected-value nodes for
#'   chi-squared or Freeman-Tukey discrepancies.
#' @param external Optional user-supplied external quantities for discrepancies
#'   that depend on inputs outside the model state.
#'
#' @return An object of class `cppp_discrepancy`.
#' @export
discrepancy <- function(kind,
                        name,
                        dataNodes = NULL,
                        modelNodes = NULL,
                        external = NULL) {
  validKinds <- c("data_only", "data_plus_model", "data_plus_external")

  if (!is.character(kind) || length(kind) != 1L) {
    stop("`kind` must be a single character string.", call. = FALSE)
  }

  if (!(kind %in% validKinds)) {
    stop(
      sprintf(
        "`kind` must be one of: %s.",
        paste(validKinds, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single character string.", call. = FALSE)
  }

  x <- list(
    kind = kind,
    name = name,
    dataNodes = dataNodes,
    modelNodes = modelNodes,
    external = external
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
#' Fill in defaults for built-in discrepancy types using information from a
#' NIMBLE model.
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

  if (out$kind == "data_only" &&
      out$name %in% c("mean", "variance")) {
    return(out)
  }

  if (out$kind == "data_plus_model" &&
      out$name %in% c("deviance")) {
    return(out)
  }

  if (out$kind == "data_plus_model" &&
      out$name %in% c("chisquared", "freemantukey")) {
    if (is.null(out$modelNodes)) {
      stop(
        sprintf(
          "Discrepancy '%s' requires `modelNodes`.",
          out$name
        ),
        call. = FALSE
      )
    }

    out$modelNodes <- model$expandNodeNames(
      out$modelNodes,
      returnScalarComponents = TRUE
    )

    if (length(out$modelNodes) != length(out$dataNodes)) {
      stop(
        sprintf(
          "`modelNodes` must match `dataNodes` in length for '%s'.",
          out$name
        ),
        call. = FALSE
      )
    }

    return(out)
  }

  stop(
    sprintf(
      paste(
        "Unsupported discrepancy combination: kind = '%s', name = '%s'."
      ),
      out$kind,
      out$name
    ),
    call. = FALSE
  )
}
