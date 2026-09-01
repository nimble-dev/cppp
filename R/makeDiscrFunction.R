#' Make offline discrepancy function
#'
#' Creates a `discFun` that computes observed and replicated discrepancies
#' for each posterior draw. This is the "offline" version used when the
#' MCMC output does not already include discrepancy columns.
#'
#' The returned function evaluates discrepancies conditional on `targetData`,
#' which may be the original observed dataset or a replicated dataset used for calibration.
#' It expects posterior samples (`MCMCSamples`) and the
#' corresponding `targetData`, and uses the model-specific functions in
#' `discConfig` to generate posterior predictive replicates and compute the
#' discrepancy values.
#'
#' @param discConfig A list containing at least two functions:
#'   \describe{
#'     \item{simulateNewDataFun}{A function
#'       `function(thetaRow, control = NULL)` that simulates one replicated
#'       dataset `y*` given a single posterior draw. The optional `control`
#'       argument can be used to pass model-specific objects (e.g., a NIMBLE
#'       model), but may be ignored if not needed.}
#'     \item{discrepancy}{A function `function(data, thetaRow, ...)` that returns
#'       the scalar discrepancy `D(data, \eqn{\theta})` for one draw.}
#'   }
#'
#' @return A function of the form
#'   `function(MCMCSamples, targetData, control = NULL, ...)` that returns a
#'   list with two numeric vectors:
#'   \itemize{
#'     \item `obs`: discrepancies computed on `targetData`.
#'     \item `sim`: discrepancies computed on posterior predictive replicates.
#'   }
#'
#' @examples
#' \dontrun{
#' discConfig <- list(
#'   simulateNewDataFun = function(thetaRow, control = NULL) rnorm(10, thetaRow),
#'   discrepancy  = function(data, thetaRow) mean((data - thetaRow)^2)
#' )
#' discFun <- makeOfflineDiscFun(discConfig)
#' }
#'
#' @export

makeOfflineDiscFun <- function(discConfig) {
  function(MCMCSamples, targetData, control = NULL, ...) {
    ## Offline discrepancy calculator: compute D(data, theta) and D(y*, theta)
    if (!is.list(discConfig) || !all(c("simulateNewDataFun", "discrepancy") %in% names(discConfig))) {
      stop("discConfig must be a list with elements 'simulateNewDataFun' and 'discrepancy'.")
    }
    simulateNewDataFun <- discConfig$simulateNewDataFun
    discrepancy  <- discConfig$discrepancy

    dObs <- apply(MCMCSamples, 1, function(thetaRow)
      discrepancy(targetData, thetaRow, ...)
    )

    dSim <- apply(MCMCSamples, 1, function(thetaRow) {
      y <- simulateNewDataFun(thetaRow, control)
      discrepancy(y, thetaRow, ...)
    })

    list(obs = dObs, sim = dSim)
  }
}
