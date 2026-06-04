#' Create a discrepancy specification
#'
#' A discrepancy specification is a lightweight user-facing description of a
#' discrepancy. The package can later complete missing defaults from a model
#' and build an R or NIMBLE implementation from the completed specification.
#'
#' @param kind Character scalar naming the structural kind of discrepancy.
#'   Current options are `"data_only"` and `"data_plus_model"`.
#' @param name Character scalar naming the discrepancy.
#' @param dataNodes Optional character vector of data node names. If `NULL`,
#'   these may be inferred from a NIMBLE model later.
#' @param modelNodes Optional character vector of model-based nodes used by the
#'   discrepancy besides the data nodes. For example, expected-value nodes for
#'   chi-squared or Freeman-Tukey discrepancies.
#'
#' @return An object of class `cppp_discrepancy`.
#' @export
discrepancy <- function(kind,
                        name,
                        dataNodes = NULL,
                        modelNodes = NULL) {
  validKinds <- c("data_only", "data_plus_model")

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
    modelNodes = modelNodes
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


#' Build an R discrepancy evaluator - this is here just for testing and building
#' but likely not needed in the future
#'
#' Create a simple R function from a completed discrepancy specification. This
#' is a lightweight way to experiment with discrepancy specifications before
#' wiring them into the NIMBLE calibration backend.
#'
#' @param model A NIMBLE model.
#' @param disc A `cppp_discrepancy` object.
#'
#' @return A function with signature `function(data = NULL)` that evaluates the
#'   discrepancy. If `data` is supplied, the evaluator temporarily inserts it
#'   into the model's data nodes, computes the discrepancy, and restores the
#'   original data values before returning.
#' @export
makeDiscrepancyRFun <- function(model, disc) {
  disc <- completeDiscrepancy(model, disc)

  function(data = NULL) {
    origData <- values(model, disc$dataNodes)
    ## on.exit tells R what to do when exiting the function
    ## I am already saying to restore original data into the model
    on.exit({
      values(model, disc$dataNodes) <- origData
      model$calculate()
    }, add = TRUE)


    if (!is.null(data)) {
      if (!is.numeric(data) || length(data) != length(disc$dataNodes)) {
        stop(
          "`data` must be a numeric vector matching the discrepancy data nodes.",
          call. = FALSE
        )
      }
      values(model, disc$dataNodes) <- data
    }

    model$calculate()

    if (disc$kind == "data_only" && disc$name == "mean") {
      return(mean(values(model, disc$dataNodes)))
    }

    if (disc$kind == "data_only" && disc$name == "variance") {
      return(stats::var(values(model, disc$dataNodes)))
    }

    if (disc$kind == "data_plus_model" && disc$name == "deviance") {
      return(2 * model$calculate(disc$dataNodes))
    }

    if (disc$kind == "data_plus_model" && disc$name == "chisquared") {
      dataVal <- values(model, disc$dataNodes)
      modelVal <- values(model, disc$modelNodes)
      return(sum((dataVal - modelVal)^2 / (modelVal + 1e-6)))
    }

    if (disc$kind == "data_plus_model" && disc$name == "freemantukey") {
      dataVal <- values(model, disc$dataNodes)
      modelVal <- values(model, disc$modelNodes)
      return(sum((sqrt(dataVal) - sqrt(modelVal))^2))
    }

    stop(
      sprintf(
        "No R evaluator implemented for kind = '%s', name = '%s'.",
        disc$kind,
        disc$name
      ),
      call. = FALSE
    )
  }
}


#' Build a NIMBLE discrepancy evaluator
#'
#' Create a NIMBLE function from a completed discrepancy specification. The
#' resulting NIMBLE function evaluates the discrepancy for the model's current
#' state.
#'
#' @param model A NIMBLE model.
#' @param disc A `cppp_discrepancy` object.
#'
#' @return A NIMBLE function object with a scalar `run()` method.
#' @export
makeDiscrepancyNimbleFun <- function(model, disc) {
  disc <- completeDiscrepancy(model, disc)

  if (disc$kind == "data_only" && disc$name == "mean") {
    return(nimbleFunction(
      setup = function(model, dataNodes) {
        modelLocal <- model
        dataNodesLocal <- dataNodes
      },
      run = function() {
        returnType(double(0))
        return(mean(values(modelLocal, dataNodesLocal)))
      }
    )(model = model, dataNodes = disc$dataNodes))
  }

  if (disc$kind == "data_only" && disc$name == "variance") {
    return(nimbleFunction(
      setup = function(model, dataNodes) {
        modelLocal <- model
        dataNodesLocal <- dataNodes
      },
      run = function() {
        returnType(double(0))
        return(var(values(modelLocal, dataNodesLocal)))
      }
    )(model = model, dataNodes = disc$dataNodes))
  }

  if (disc$kind == "data_plus_model" && disc$name == "deviance") {
    return(nimbleFunction(
      setup = function(model, dataNodes) {
        modelLocal <- model
        dataNodesLocal <- dataNodes
      },
      run = function() {
        returnType(double(0))
        return(2 * modelLocal$calculate(dataNodesLocal))
      }
    )(model = model, dataNodes = disc$dataNodes))
  }

  if (disc$kind == "data_plus_model" && disc$name == "chisquared") {
    return(nimbleFunction(
      setup = function(model, dataNodes, modelNodes) {
        modelLocal <- model
        dataNodesLocal <- dataNodes
        modelNodesLocal <- modelNodes
      },
      run = function() {
        dataVal <- values(modelLocal, dataNodesLocal)
        modelVal <- values(modelLocal, modelNodesLocal)
        returnType(double(0))
        return(sum((dataVal - modelVal)^2 / (modelVal + 1e-6)))
      }
    )(model = model, dataNodes = disc$dataNodes, modelNodes = disc$modelNodes))
  }

  if (disc$kind == "data_plus_model" && disc$name == "freemantukey") {
    return(nimbleFunction(
      setup = function(model, dataNodes, modelNodes) {
        modelLocal <- model
        dataNodesLocal <- dataNodes
        modelNodesLocal <- modelNodes
      },
      run = function() {
        dataVal <- values(modelLocal, dataNodesLocal)
        modelVal <- values(modelLocal, modelNodesLocal)
        returnType(double(0))
        return(sum((sqrt(dataVal) - sqrt(modelVal))^2))
      }
    )(model = model, dataNodes = disc$dataNodes, modelNodes = disc$modelNodes))
  }

  stop(
    sprintf(
      "No NIMBLE evaluator implemented for kind = '%s', name = '%s'.",
      disc$kind,
      disc$name
    ),
    call. = FALSE
  )
}
