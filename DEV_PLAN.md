# Development plan

Last updated: 2 October 2026.

## How the package works

- The user writes one or more `discrepancy()` and one `simulation()`.
- The package builds a **calculator**. For each posterior draw it gives two
  numbers per discrepancy: one for the data, one for a replicate.
- `runCalibration()` turns those numbers into a PPP and then a CPPP.

## What works

- Five built-in discrepancies. A user discrepancy is a nimbleFunction with
  `contains = discrepancyBase`, a setup `(model, dataNodes, modelNodes)` and a
  `run()` returning one number.
- `simulation(mode, latentNodes)`. Conditional redraws the data. Marginal
  redraws the named latents and everything below them.
- `completeNodes()` works out every node list once. The calculator, the
  simulate function and the MCMC all use it.
- `runCalibrationNIMBLE()` compiles the model, the MCMC and the calculator in
  one call, on one model.
- The calculator saves and restores only the nodes it changes.
- `makeDiscrepancyExtractor()` reads discrepancies the MCMC already computed.
  Placeholder: nothing feeds it yet.

## Done recently

- `78a01ef` `completeNodes()` and `simulation(mode, latentNodes)`.
- `ee625bd` Calculator, simulate function and `runCalibrationNIMBLE()` use
  `completeNodes()`. The calculator restores the model with `modelValues`.
- `afe7b6c` Default `control$disc` entries filled from the nodes.

## Decided, not yet done

- **Three layers.**
  - `runCalibration()`: the engine, for models written in plain R.
  - `runCalibrationNIMBLE()`: the advanced NIMBLE function.
  - `cppp_nimble()`: a simple wrapper for new users. Later.
- **No user simulation function in `runCalibrationNIMBLE()`.** Conditional and
  marginal cover the real cases. Drop `simulateNewDataFun`. Anything else goes
  through `runCalibration()`.
- **R discrepancies use the package's replicates.** The package builds the
  simulate function from `simulation()` and passes it to the user's `discFun`.
- **Same replicate for every discrepancy.** For each draw, all discrepancies
  see the same replicate dataset.
- **`print()`, `summary()`, `plot()`** for the result.
- **Update the docs** (README, cheatsheet) for `simulation(mode, latentNodes)`.

## Open questions

- **One `discrepancies` argument for both kinds?** nimbleFunction and R
  function, told apart by the type of `fun`. `discFun` would then go too. What
  does an R discrepancy receive: `function(data, theta)`? Waiting for input.
- **The `control` list.** It carries duplicates (names and nodes), an
  uncompiled model, and silently ignores top-level entries when sublists
  exist. Its shape depends on the question above.
- **Prior predictive.** If wanted, a third mode of `simulation()`.

## Next jobs

1. **Keep the shape of the nodes.** Started: only data, parameters, latents
   and the discrepancy nodes are still flattened, all in `simulationSpec.R`
   and `discrepancySpec.R`. Data and parameters need it. Check the other
   three.
2. **Run calibration replicates in parallel.** One set of compiled objects per
   core, built once. `future` is the suggested tool. The loop is in
   `runCalibration.R`.
3. **Reuse the long chain's sampler tuning** in the short chains. Daniel Turek
   has code on the nimble forum. Not the same as `transferAutocorrelation()`.
4. **An end-to-end test in marginal mode.**

## Rules worth remembering

- After setting a draw, calculate the parameters' dependencies, never the
  parameters (`self = FALSE`). A parameter can be derived (`sigma` from
  `log_sigma`), and calculating it would overwrite the drawn value.
- A discrepancy's data nodes must be among the data nodes.
- Compile everything in one `compileNimble(list(...))` call.
- Draws and replicates travel between pieces as R matrices and vectors. Using
  `modelValues` there was considered and dropped: the copying is tiny next to
  the MCMC, and the engine works with R matrices. `modelValues` is used only
  inside the calculator, to save and restore the model.
- *Calculator* computes discrepancies after the MCMC. *Extractor* reads values
  the MCMC computed. Neither makes the PPP: that is `runCalibration()`.
- Not everything compiles (`sort()` does not). Wrap R code with
  `nimbleRcall()`. One defined inside the package must be exported.
- Compiled code cannot pick columns by name. The wrapper orders the parameter
  columns first.
- `asymm` in the Newcomb example hardcodes order statistics 6 and 61 (n = 66).

## Files

| File | Holds |
|---|---|
| `R/discrepancySpec.R` | `discrepancy()`, `completeDiscrepancy()`, `makeDiscrepancyNimbleFun()` |
| `R/builtinDiscrepancies.R` | `discrepancyBase` and the five built-ins |
| `R/simulationSpec.R` | `simulation()`, `completeNodes()` |
| `R/makeDiscrepancyCalculator.R` | the calculator |
| `R/makeSimulateNewDataFun.R` | makes one replicate from a draw |
| `R/makeDiscrepancyExtractor.R` | reads discrepancies from MCMC output |
| `R/makeDiscrFunction.R` | `makeOfflineDiscFun()`, the plain-R route |
| `R/runCalibration.R` | the engine |
| `R/runCalibrationNIMBLE.R` | the NIMBLE wrapper |
| `R/cpppResult-class.R` | the result object |
| `R/transferAutocorrelation.R` | placeholder |
| `inst/examples/newcomb_spec_offline.R` | worked example |
