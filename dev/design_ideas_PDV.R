# Sketch of ideas

runCalibration <- function(MCMC_samples, # assume the main run happened already
                           MCMC_fun,
                           new_data_fun,

                           disc_fun

                           ) {
  ## No nimble-specific concepts or code in this function.
  #
  ## Ignore where the "main" happens
  ppp <- list()
  for(i in 1:num_reps) {
    new_data <- new_data_fun(MCMC_samples)
    new_samples <- MCMC_fun(new_data, control)
    discrepancies <- disc_fun(new_samples, ...)
    # ...
  }
}

# Two options for what disc_fun would do:
# If the discrepancies were really computed onine during MCMC_fun,
# they may just be in an output column of new_samples.
# Then disc_fun could simply be a function that picks out the right column.
#
#

make_col_disc_fun <- function(col_name) {
  function(MCMC_samples) {
    MCMC_samples[, col_name]
  }
}

make_offline_disc_fun <- function(control) {
  function(MCMC_samples) {
    # Actually runs posterior predictive simulations etc.
  }
}

make_MCMCfun <- function(niter, data_nodes, cmodel, cmcmc) {
  function(new_data, control) {
    cmodel[[data_nodes]] <- new_data
    samples <- nimble::runMCMC(cmcmc,niter=niter )
  }

}

# It's fully general, because we'll just call "disc_fun" with some sufficiently general set of arguents (MCMC_samples)
# Any other problem-specific controls can get baked in my the "make_...fun" step (e.g. col_name, or control)


runCalibrationNIMBLE <- function(nimble_specific_args) {
  # Use all the nimble-specific steps to create the generic inputs
  #  that are a set of functions for runCalibration
  # build and compile the model and mcmc pieces
  # Set up "MCMC_fun"
  MCMCfun <- make_MCMCfun(calib_niter, simulated_data_nodes, cmodel, cmcmc)
  # Set up disc_fun
  # Set up the data_sim_fun
  disc_fun <- make_col_disc_fun("some_disc_output")
  results <- runCalibration(MCMCfun, disc_fun)
  results

}

# Maybe we or someone else will come along later and create "runCalibration_other"

# To think about: how to make it easy for a user to provide a discrepancy function with a set of nodes to simulate and plug it in easily without diving into lots of internal details.

#Example

new_disc <- function(y, mu) {
  return(sum((y-mu)^2))
}

# But maybe we will defer this to avoid getting to tricky right now!
