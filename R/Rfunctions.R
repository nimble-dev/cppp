# R/Rfunctions.R

#' Compute the calibrated posterior predictive p-value
#'
#' @param obsPPP Numeric scalar in \eqn{[0,1]}. The observed PPP.
#' @param repPPP Numeric vector in \eqn{[0,1]}. PPPs from calibration datasets.
#'
#' @return Numeric scalar in \eqn{[0,1]}, the calibrated PPP.
#' @export
computeCppp <- function(obsPPP, repPPP) {
  # basic validation
  if (!is.numeric(obsPPP) || length(obsPPP) != 1L)
    stop("obsPPP must be a numeric scalar.", call. = FALSE)
  if (!is.numeric(repPPP))
    stop("repPPP must be numeric.", call. = FALSE)

  # handle empty vector as NA ??
  if (length(repPPP) == 0L)
    return(NA_real_)

  mean(repPPP <= obsPPP)
}


#' Estimate transfer autocorrelation (placeholder)
#'
#' Placeholder for the function that will estimate the transfer
#' autocorrelation from the observed Delta-chain
#'
#' @param deltaChain Data or structure containing \eqn{D(y^{rep}, \theta) - D(y^{obs}, \theta)} values.
#' @param obsPPP Numeric scalar in  \eqn{[0,1]}.
#' @param repPPP Numeric vector in  \eqn{[0,1]}.
#' @param mTilde Numeric vector or scalar; effective Monte Carlo sample size per replication.
#' @param ... Additional arguments (future use).
#'
#' @return List with components to be defined (e.g., `iact`, `essTransfer`).
#' @export
transferAutocorrelation <- function(deltaChain, obsPPP, repPPP, mTilde, ...) {
  stop("Not implemented", call. = FALSE)
}
