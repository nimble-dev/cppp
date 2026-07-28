# cppp

A package for Calibrated posterior predictive p-values.

> ### ⚠️ Work in progress — not ready for use
>
> This is an active research project in early development. It is public so
> collaborators can follow along, not because it is ready.
>
> - The interface changes often, without notice or deprecation.
> - There is no released version, and nothing here should be depended on.
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

### Main functions 
| Function | Purpose |
|-----------|----------|
| `runCalibration()` | Generic function for calibration |
| `runCalibrationNIMBLE()` | NIMBLE-specific setup that builds the generic inputs and calls `runCalibration()`. |

### Helpers
| Function | Role |
|-----------|------|
| `discrepancy()` | Describes one discrepancy: its name, the nodes it reads, and optionally your own nimbleFunction. |
| `simulation()` | Describes how replicated datasets are generated. |
| `makeColDiscFun()` | Returns a function that extracts an “online” discrepancy column. |
| `makeOfflineDiscFun()` | Returns a function that computes discrepancies “offline.” |
| `computeCppp()` | Computes cppp from $\hat p_{\text{obs}}$ and $\hat p_j$. |
| `transferAutocorrelation()` | Transfer-ESS machinery for the cppp variance. Placeholder — not implemented. |

### Data object

`cpppResult` (S3), currently holding:

  - `CPPP` — one value per discrepancy
  - `obsPPP`, `repPPP` — observed and replicated posterior predictive p-values
  - `discrepancies` — the observed and replicated discrepancy values
  - `drawnIndices` — which posterior draws seeded the calibration replicates

Standard errors, confidence intervals, and `print` / `summary` / `plot` methods
are planned but not implemented.

---

## Typical workflow 

1. **Build your MCMC engine**
   - For NIMBLE: configure a model and compiled MCMC object.
   - Otherwise: supply a function `MCMCFun(targetData, control)` that runs an MCMC and returns samples.

2. **Provide a data simulator**
   - `simulateNewDataFun(thetaRow, control)` draws a dataset $\tilde{y}$ from the model posterior predictive.

3. **Describe your discrepancies**
   - `discrepancy("mean")` names one the package ships; `discrepancy("asymm", modelNodes = "mu", fun = myAsymm)` supplies your own as a nimbleFunction.
   - *Online*: the discrepancy or PPP is computed during MCMC; just extract the column.
   - *Offline*: compute $D(y,\theta)$ and $D(y^*,\theta)$ after sampling.
   - Note: the step that turns these descriptions into what `runCalibration()` consumes is still being built.

4. **Run calibration**
   - Call `runCalibration()` with your functions and number of replicates $r$.
   - the function iteratively: simulate → run short chain → compute $\hat{p}_j$.

5. **Compute cppp and variance**
   - `computeCppp(obsPPP, repPPP)`
   - `transferAutocorrelation(deltaChain, obsPPP, repPPP, mTilde)` — not implemented yet

---

## Example

TO DO 
