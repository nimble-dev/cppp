#' Read discrepancies the MCMC already computed
#'
#' Placeholder, until the online route settles. When the discrepancies are
#' computed during the MCMC run as a derived quantity, the values are already
#' in the output and this reads them out, in the shape
#' [makeDiscrepancyCalculator()] returns.
#'
#' With the stem `"discrepancy"` the columns read are `discrepancy_model` (the
#' observed side) and `discrepancy_simulated` (the replicate side). The naming
#' is not fixed yet, which is why the stem and the suffixes are arguments.
#'
#' @param discrepancies Character vector of column stems. Each stem also names
#'   the matching output column.
#' @param obsSuffix,simSuffix Character, appended to each stem to give the
#'   observed and the replicated column.
#'
#' @return A function `function(MCMCSamples, ...)` returning `list(obs, sim)`,
#'   each a draws x discrepancies matrix.
#' @export
makeDiscrepancyExtractor <- function(discrepancies = "discrepancy",
                                     obsSuffix = "_model",
                                     simSuffix = "_simulated") {

  obsCols <- paste0(discrepancies, obsSuffix)
  simCols <- paste0(discrepancies, simSuffix)

  function(MCMCSamples, ...) {
    MCMCSamples <- as.matrix(MCMCSamples)

    absent <- setdiff(c(obsCols, simCols), colnames(MCMCSamples))
    if (length(absent) > 0L) {
      stop("`MCMCSamples` has no column for: ",
           paste(absent, collapse = ", "), call. = FALSE)
    }

    obs <- MCMCSamples[, obsCols, drop = FALSE]
    sim <- MCMCSamples[, simCols, drop = FALSE]
    colnames(obs) <- colnames(sim) <- discrepancies

    list(obs = obs, sim = sim)
  }
}
