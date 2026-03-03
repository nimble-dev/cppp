#' Create a cpppResult object
#'
#' Constructs an S3 object of class `cpppResult`, which stores the results
#' related to the calibration of the predictive p-value.
#'
#' @param CPPP Numeric scalar. The approximated calibrated posterior predictive p-value.
#' @param repPPP Numeric vector. The posterior predictive p-values across calibration datasets.
#' @param obsPPP Numeric vector. The observed posterior predictive p-values for the observed data.
#'
#' @param discrepancies Optional list or array containing discrepancy values for
#'   each replication and dataset (observed and simulated).
#'
#' @param ... Additional elements to store in the object.
#' SP E.g., do we want to save in this object functions used in the pipeline?
#'
#' @return An object of class `cpppResult`.
#' @export
newCpppResult <- function(CPPP = NA_real_,
                           repPPP = numeric(),
                           obsPPP = numeric(),
                           discrepancies = NULL,
                           ...) {
  x <- list(
    CPPP = CPPP,
    repPPP = repPPP,
    obsPPP = obsPPP,
    discrepancies = discrepancies,
    ...
  )
  class(x) <- c("cpppResult", "list")
  validateCpppResult(x)
}

# Minimal internal validator: keeps inputs sane but stays permissive

#' @keywords internal
validateCpppResult <- function(x) {
  stopifnot(inherits(x, "cpppResult"))

  # CPPP: allow NA for now; if finite, must be a scalar in [0,1]
  if (!(length(x$CPPP) == 1L && is.numeric(x$CPPP))) {
    stop("`CPPP` must be a numeric scalar (NA allowed).", call. = FALSE)
  }
  if (is.finite(x$CPPP) && (x$CPPP < 0 || x$CPPP > 1)) {
    stop("`CPPP` must be in [0,1] when finite.", call. = FALSE)
  }

  # repPPP: numeric vector, all finite in [0,1] if provided
  if (!is.numeric(x$repPPP)) {
    stop("`repPPP` must be a numeric vector.", call. = FALSE)
  }
  if (length(x$repPPP) && (!all(is.finite(x$repPPP)) || any(x$repPPP < 0 | x$repPPP > 1))) {
    stop("All `repPPP` values must be finite and in [0,1].", call. = FALSE)
  }

  # obsPPP: numeric vector, all finite in [0,1] if provided
  if (!is.numeric(x$obsPPP)) {
    stop("`obsPPP` must be a numeric vector.", call. = FALSE)
  }
  if (length(x$obsPPP) && (!all(is.finite(x$obsPPP)) || any(x$obsPPP < 0 | x$obsPPP > 1))) {
    stop("All `obsPPP` values must be finite and in [0,1].", call. = FALSE)
  }

  # discrepancies: optional; no strict contract yet
  # (kept permissive to allow list/array/tibble as you iterate)

  invisible(x)
}
