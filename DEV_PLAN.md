# Where the discrepancy work stands

Last worked on: 1 September 2026. Everything below is current.

---

## The picture in one paragraph

You write **specifications**: one or more `discrepancy()` and one
`simulation()`. The package turns those into a **calculator**, which walks
through the posterior draws and produces two numbers per draw per discrepancy —
one for your dataset, one for a replicate simulated from that draw. Those go to
`runCalibration()`, which turns them into a PPP and then a CPPP. The
calculation of the discrepancies and the calculation of the PPP are two
separate jobs, and the names say so.

---

## What works today

- **Specifications.** `discrepancy(name, dataNodes, modelNodes, fun)` and
  `simulation(mode, dataNodes, simulateNodes)`. Five built-in discrepancies
  ship with the package; naming one is a shortcut for supplying it yourself.
- **The shape of a discrepancy.** A nimbleFunction with
  `contains = discrepancyBase`, a setup taking `(model, dataNodes, modelNodes)`
  and a `run()` with no arguments returning one number. Built-ins and your own
  are the same thing.
- **`makeDiscrepancyCalculator()`.** The main piece. The loop over draws lives
  inside a nimbleFunction, with the discrepancies held in a
  `nimbleFunctionList`, so the whole thing compiles in one call and the loop
  runs in compiled code. Returns a function of
  `(MCMCSamples, targetData, control)` giving `list(obs, sim)`, each one
  draws x discrepancies with the discrepancy names on the columns.
- **The calculator is wired into `runCalibrationNIMBLE()`.** Pass
  `discrepancies =` and `simulation =` and it builds the calculator and the
  replicate-simulating function for you. The older `discFun` argument still
  works.
- **`makeDiscrepancyExtractor()`** reads several named discrepancy columns
  straight out of MCMC output, in the same shape the calculator returns.
- **The engine.** `runCalibration()` and `runCalibrationNIMBLE()` handle
  several discrepancies: `obsPPP` and `CPPP` come back as named vectors,
  `repPPP` as a matrix.

Tests: 97 passing. `devtools::test()`.

---

## What to do next

Four jobs. The first two are groundwork, the third needs them, the fourth is
separate.

1. **Keep the shape of the nodes.** Every place that handles node names
   flattens them with `expandNodeNames(..., returnScalarComponents = TRUE)`, so
   a vector `y` of length 66 becomes 66 loose strings. The shape is what tells
   us how long an output should be and how to read values back out in the right
   arrangement, so we should keep it. Decide what to hold beside the flat names,
   and which places need the flat names at all — writing values into the model
   does, dependency lookups and monitors probably do not. Five files do it, and
   the note at `makeDiscrepancyCalculator.R` line 64 already wonders about one.

2. **One model, not three.** A call to `runCalibrationNIMBLE()` makes three
   copies and compiles twice: one for the calculator, one for the simulate
   function — never compiled, so it runs in R once per replicate — and the
   original for the MCMC. They never run at the same time, so one is enough.
   Compile everything against it in a single `compileNimble(list(...))`, which
   should also drop `resetFunctions = TRUE`, and move data between the pieces
   with `modelValues` rather than R matrices. Watch the starting values: the
   calculator leaves the model holding the last draw it looked at.

3. **Run the calibration replicates in parallel.** They go one after another in
   a plain loop, `runCalibration.R` line 141. Build a small collection of
   compiled nimble objects — model, MCMC, calculator, simulate function — one
   per core, and let each core keep its own across one replicate after another.
   With 100 replicates on 5 cores that is 5 compilations, not 100. Compiled
   objects cannot be sent between processes, so each core builds its own; the
   point is that it builds it once. `future` is the suggested tool. Needs job 2.

4. **Let the short chains reuse the long chain's sampler tuning.** Every short
   replicate chain re-tunes step sizes and proposal covariances from scratch,
   eating a burn-in that is already short. Start them from what the long chain
   settled on instead. Daniel Turek has code on the nimble forum for reading
   those values out and putting them back — find it first, it decides the rest.
   Not the same as `transferAutocorrelation()`, which is about the delta chain.

---

## Decisions, and why

**One model, shared.** Earlier this plan said the calculator needed its own
copy because the MCMC is constantly changing the model. Reversed: the MCMC and
the calculator never run at the same time, so there is nothing to collide.
The calculator writes each draw in before it reads anything, so it does not
care what the MCMC left behind. See task 2 for the one place this shows up.

**Marginal simulation: you name the nodes.** `simulation("marginal")` needs
`simulateNodes` — the latent nodes plus the data nodes. There is no default.
We built automatic derivation (everything below the parameters) and took it
out: which nodes count as latent depends on what you chose to treat as a
parameter, and `runCalibrationNIMBLE`'s own default sweeps latent states in
with the parameters, so marginal quietly did conditional's job. Worth
revisiting only if we find a reliable way to tell the two apart.

**Calculate below the parameters, never the parameters themselves.** After
setting a draw, the calculator calls
`model$calculate(getDependencies(paramNodes, self = FALSE))`. This matters
because a parameter can be a *derived* node — `sigma`, when the prior is on
`log(sigma)`. Calculating a derived node rewrites its value from its parents,
so a whole-model `calculate()` would throw away the drawn `sigma` and freeze it
across every draw, silently. With `self = FALSE`, monitoring `sigma` and
monitoring `log_sigma` both work. There is a test for exactly this.

**A discrepancy's `dataNodes` must sit inside the simulation's.** The dataset
is written into the simulation's nodes; each discrepancy reads its own. If a
discrepancy reads something the simulation never writes, it is reading
leftovers. The calculator checks and stops. A discrepancy reading only part of
the data is still allowed.

**Compile in one call.** NIMBLE only handles adding nimbleFunctions to an
existing project in limited circumstances — it reuses compiled pieces and gives
up if they lack what it needs. The manual's own advice is to compile several
units together, which is what the `nimbleFunctionList` buys us:
`compileNimble(list(model, calcNF))`. `resetFunctions = TRUE` is the other way
out, and is why it currently appears in `runCalibrationNIMBLE`; task 2 should
remove the need for it.

**Names.** *Calculator* computes discrepancy values from the model, after the
MCMC. *Extractor* reads values the MCMC already computed. Neither turns them
into a PPP — that is `runCalibration()`'s job, and keeping the names apart keeps
that boundary visible.

---

## Gotchas

- **`contains = discrepancyBase` is required** on every discrepancy. Without
  it, it cannot go into the list.
- **Not everything compiles.** `sort()` does not. Wrap R functions with
  `nimbleRcall()`, as `sortR` does in the Newcomb example. A `nimbleRcall`
  defined inside the package must be exported or NIMBLE will not find it; ones
  you define in your own script are fine.
- **Give the calculator an uncompiled model.** It compiles its own copy. A
  compiled model keeps a link back at `cmodel$Rmodel` if you ever need it.
  (Task 2 changes this.)
- **Compiled code cannot pick out columns by name.** The wrapper picks out the
  parameter columns in `paramNodes` order before handing the matrix over.
- **`asymm` in the Newcomb example hardcodes order statistics 6 and 61**, which
  are only right for n = 66.
- Run `devtools::document()` after adding anything exported.

---

## Files

| file | what it holds |
|---|---|
| `R/discrepancySpec.R` | `discrepancy()`, `completeDiscrepancy()`, `makeDiscrepancyNimbleFun()` |
| `R/builtinDiscrepancies.R` | `discrepancyBase` and the five built-ins |
| `R/simulationSpec.R` | `simulation()`, `completeSimulation()` |
| `R/makeDiscrepancyCalculator.R` | the calculator and the nimbleFunction holding the loop |
| `R/makeDiscrepancyExtractor.R` | reads discrepancies the MCMC already computed |
| `R/makeSimulateNewDataFun.R` | simulates one replicate dataset from a draw |
| `R/makeDiscrFunction.R` | `makeOfflineDiscFun`, the older plain-R route kept for cross-checking |
| `R/runCalibration.R` | the generic engine — no model, just draws and functions |
| `R/runCalibrationNIMBLE.R` | the NIMBLE wrapper |
| `R/cpppResult-class.R` | the result object |
| `R/transferAutocorrelation.R` | placeholder, not implemented |
| `inst/examples/newcomb_spec_offline.R` | the worked example |
