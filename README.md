# cppp

A package for Calibrated posterior predictive p-values.

> ### ⚠️ Work in progress — not ready for use
>
> This is an active research project in early development. It is public so
> collaborators can follow along, not because it is ready.
>
> There is no released version, and nothing here should be depended on.
>
> If you want to try it, get in touch first.

It provides a general framework for a MCMC engine—to:

1. Compute calibrated posterior p-values (cppp),
2. Estimate their Monte Carlo variance using the idea of the **transfer effective sample size (ESS)**
3. Can handle different MCMC engines. **NIMBLE** and R for now, other MCMC engines later.

---

## Concept

Given data $y$, a model $p(\theta \mid y) \propto p(y \mid \theta) \pi (\theta)$, and a discrepancy function $D(y,\theta)$:

1. Run an long MCMC chain to obtain draws from the posterior $p(\theta \mid y)$. With $M$ draws, we sample new datafrom the posterior predictive of the data $p(y^* \mid \theta_i)$ and compute 

$$
 \Delta_i = D(y^*_i, \theta_i) - D(y, \theta_i),
$$

This chain of $\Delta = \{ \Delta_i \}$ collects the observed discrepancies.

2. Generate $r$ *calibration replicates*:
   - simulate new datasets $\tilde{y}_j$ from the model,
   - run short chains of length $\tilde{m} $ for $p(\theta \mid \tilde{y}_j)$,
   - compute short-run posterior-predictive p-values $\hat{p}_j$.

3. Estimate the CPPP:

$$
 \widehat{\text{cppp}}
 = \frac{1}{r}\sum_{j=1}^r
   \mathbf{1}\{\hat{p}_j \le \hat{p}_{\text{obs}}\}.
$$

4. Estimate the Monte Carlo variance using the *transfer ESS* idea:
   match each $\hat{p}_j$ to a quantile on the observed  $\Delta$ -chain,
   compute the transfer autocorrelation, and estimate the cppp variance.
  
---

## Package architecture

There is a one-page architecture map: the two routes through the package, what
you write, what comes back, and what is not not working yet.

**[View the architecture map](https://raw.githack.com/nimble-dev/cppp/master/inst/cheatsheet/cppp-architecture.html)**
&nbsp;·&nbsp;
[source](inst/cheatsheet/cppp-architecture.html)

Or, once the package is installed:

```r
browseURL(system.file("cheatsheet", "cppp-architecture.html", package = "cppp"))
```

### Main functions

| Function | Purpose |
|-----------|----------|
| `runCalibration()` | The generic engine. No model in it — just posterior draws and the functions you give it. |
| `runCalibrationNIMBLE()` | The NIMBLE wrapper. Give it a model and your specifications and it builds the rest itself. |

### Describing what you want

| Function | Role |
|-----------|------|
| `discrepancy()` | Describes one discrepancy: this is a list containing name, the nodes it reads, and optionally a custom `nimbleFunction`. |
| `simulation()` | Describes how a replicate dataset is generated. Two options: `"conditional"` (redraw the data) or `"marginal"` (redraw named latent nodes too). |
| `discrepancyBase` | a virtual base class. A custom nimbleFunction discrepancy must `contain` it. |

Five discrepancies are builtin with the package, and naming one is a shortcut for
writing it yourself: `mean`, `variance`, `deviance`, `chisquared`,
`freemantukey`.

#### Implementation notes

Both simulation() and discrepancy() take a dataNodes argument. The simulation's set has to cover the discrepancy's:

```
sim <- simulation("conditional", dataNodes = c("y", "z"))  # writes y and z

discrepancy("mean", dataNodes = "y")           # reads y  -- fine, a subset
discrepancy("mean", dataNodes = c("y", "z"))   # reads both -- fine
discrepancy("mean", dataNodes = "w")           # reads w  -- error
```

The simulation writes first and the discrepancy reads afterwards. A node the simulation never writes still holds the value from the previous draw, so a discrepancy reading it would not actually use a new replicate. The calculator checks the two sets before it starts and stops, naming the nodes at fault.

For `paramNodes` the story is the same. Before either one does anything, the calculator and the simulator each write a posterior draw into their model. Both look up the values they need by node name, so both must be given the same parameter nodes — and every one of those nodes has to appear as a column in the posterior samples. If one is missing, you get an error naming it.

`runCalibrationNIMBLE()` arranges all of this for you. It passes your dataNames to the simulation and your paramNames to both the calculator and the simulator, so they cannot disagree. You only need to think about it if you build a calculator or a simulator yourself.

**Each piece gets its own copy of the model.** Three copies are made in all:

| copy | used by | compiled? |
|---|---|---|
| the original | the MCMC | yes |
| `model$newModel()` | the discrepancy calculator | yes, together with the discrepancies |
| `model$newModel()` | the replicate simulator | no, it stays plain R |

The three run one after another, never at the same time: simulate a replicate
dataset, run a short chain on it, then compute the discrepancies. What the
copies avoid is not a clash but the state each one leaves behind. A model
remembers the last values written into it, and the MCMC starts its next chain
from wherever its model happens to be. Separate copies keep one step from
moving another's starting point.

The cost worth knowing about is that the model is compiled twice, once for the
MCMC and once with the discrepancies. That is the slow part of setting up a
run. The simulator's copy is cheap by comparison, since it is never compiled.

> **Open question.** The third copy may not be needed. Once the calculator is
> compiled it works on the C model and never touches the R model it was built
> from, so the simulator could reuse that one. The saving is small and the risk
> is a state bug that gives wrong numbers rather than an error, so this is
> written down for discussion rather than acted on.

**Watch out for derived parameters.** Say the prior is on `log(sigma)` but you
monitor `sigma`. After writing a draw into the model, the calculator
recalculates only the nodes *below* the parameters. It never recalculates the
parameters themselves. If it did, it would rebuild `sigma` from its parents and
throw away the value just drawn — the same `sigma` for every draw, with no error
to warn you. Because it skips them, monitoring `sigma` and monitoring
`log_sigma` both work.

### Turning descriptions into numbers

| Function | Role |
|-----------|------|
| `makeDiscrepancyCalculator()` | Builds the calculator that simulate from the posterior predictive and returns $D(y,\theta)$ and $D(y^*,\theta)$ per each replicated dataset, per discrepancy. Discrepancies and loop over the discrepancies compile together in one call using `makeDiscrepancyNimbleFun()`. |
| `makeDiscrepancyNimbleFun()` | Turns one specification into a nimbleFunction. Mostly used by the calculator. |
| `makeSimulateNewDataFun()` | Builds the function that simulates one replicate dataset from a draw. |
| `makeDiscrepancyExtractor()` | **Placeholder.** The online counterpart — reads discrepancies the MCMC already computed instead of recomputing them. |

`runCalibrationNIMBLE()` calls the first two for you when you pass
`discrepancies` and `simulation`, so most of the time you never touch them
directly.

### Results

`runCalibration()` returns a list with:

  - `CPPP` — one value per discrepancy
  - `obsPPP`, `repPPP` — observed and replicated posterior predictive p-values
  - `discrepancies` — the observed and replicated discrepancy values
  - `drawnIndices` — which posterior draws seeded the calibration replicates

`runCalibration()` returns this as a `cpppResult` S3 object, built and checked
by `newCpppResult()`. Standard errors, confidence intervals, and
`print` / `summary` / `plot` methods are planned but not implemented, so for
now it prints as a plain list.

### Placeholders

Two pieces are deliberately unfinished, and neither is tested:

| File | What it is |
|---|---|
| `R/transferAutocorrelation.R` | `transferAutocorrelation()` stops with "Not implemented". Corresponds to step 4 of the concept above, i.e.,  the transfer-ESS estimate of the cppp variance. |
| `R/makeDiscrepancyExtractor.R` | `makeDiscrepancyExtractor()` reads discrepancy columns that NIMBLE computed during the MCMC run. It is written against the current derived-quantity output (`discrepancy_model`, `discrepancy_simulated`), but that format is not settled, and nothing is wired up to feed it yet. |

---

## Typical workflow

1. **Build the model** — an ordinary NIMBLE model with your data in it.

2. **Describe your discrepancies** — `discrepancy("mean")` names one the package
   ships; `discrepancy("asymm", modelNodes = "mu", fun = myAsymm)` supplies your
   own as a nimbleFunction with `contains = discrepancyBase`.

3. **Describe the replicates** — `simulation("conditional")` redraws the data
   from each posterior draw. `simulation("marginal", simulateNodes = ...)`
   redraws latent nodes as well, and you must name them.

4. **Run the calibration** — hand both to `runCalibrationNIMBLE()`, which runs
   the long chain, simulates $r$ replicate datasets, runs a short chain on each,
   and turns the results into a PPP and a CPPP.

5. **Estimate the variance** — `transferAutocorrelation()`, once implemented.

---

## Example

The Newcomb light-speed measurements, with two discrepancies at once: `mean`,
which a normal model should fit, and an asymmetry statistic, which it should
not. The full script is in
[`inst/examples/newcomb_spec_offline.R`](inst/examples/newcomb_spec_offline.R).

```r
library(nimble)
library(cppp)

lightPath   <- system.file("examples", "light.txt", package = "cppp")
newcombData <- list(y = read.table(lightPath)$V1)

newcombModel <- nimbleModel(
  code = nimbleCode({
    for (i in 1:n) y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dflat()
    log(sigma) ~ dflat()
  }),
  data      = newcombData,
  constants = list(n = length(newcombData$y)),
  inits     = list(mu = 0, log_sigma = 2)
)

## Your own discrepancy: how lopsided the two tails are around mu.
## sort() does not compile, so wrap R's.
sortR <- nimbleRcall(prototype = function(x = double(1)) {},
                     Rfun = "sort", returnType = double(1))

asymmDisc <- nimbleFunction(
  contains = discrepancyBase,
  setup = function(model, dataNodes, modelNodes) {},
  run = function() {
    ys <- sortR(values(model, dataNodes))
    mu <- values(model, modelNodes)[1]
    returnType(double(0))
    return(abs(ys[61] - mu) - abs(ys[6] - mu))
  }
)

res <- runCalibrationNIMBLE(
  model         = newcombModel,
  dataNames     = "y",
  paramNames    = c("mu", "sigma"),
  discrepancies = list(discrepancy("mean"),
                       discrepancy("asymm", modelNodes = "mu", fun = asymmDisc)),
  simulation    = simulation("conditional"),
  nReps         = 20
)

res$obsPPP   # one per discrepancy; mean near 0.5, asymm far from it
res$CPPP     # one per discrepancy
res$repPPP   # nReps x 2
```
