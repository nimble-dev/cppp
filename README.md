# cpppPrototype

This package implements the calibration procedure described in  
*Computational Methods for Fast Bayesian Model Assessment via Calibrated Posterior p-values*.

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

3. Combine:

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
| `make_col_disc_fun()` | Returns a function that extracts an “online” discrepancy column. |
| `make_offline_disc_fun()` | Returns a function that computes discrepancies “offline.” |
| `compute_cppp()` | Computes cppp from $\hat p_{\text{obs}}$ and $\hat p_j$. |
| `transfer_ess_variance()` | Estimates cppp variance via transfer ESS. |

Other possible helpers? 

| Function | Role |
|-----------|------|
| `make_data_sim_fun()` | Returns a function that simulates new datasets $\tilde y$. |
| `make_MCMCfun()` | [need to check what that was] |

### Data object

`cpppResult` (S3) 

  - cppp estimate, se, confidence interval
  - observed ppp, replicated ppp
  - observed discrrepancies, replicated discrepancies (optional?)
  - informations about the calibration procedure? number of replicates and mcmc per replicated
  - methods: `print`, `summary`, `plot`

---

## Typical workflow 

1. **Build your MCMC engine**
   - For NIMBLE: configure a model and compiled MCMC object.
   - Otherwise: supply a function `MCMC_fun(new_data, control)` that runs an MCMC and returns samples.

2. **Provide a data simulator**
   - `new_data_fun()` draws a dataset $\tilde{y}$ from the model posterior predictive.

3. **Provide a discrepancy handler**
   - *Online*: the discrepancy or PPP is computed during MCMC; just extract the column.
   - *Offline*: compute $D(y,\theta)$ and $D(y^*,\theta)$ after sampling.

4. **Run calibration**
   - Call `runCalibration()` with your functions and number of replicates $r$.
   - the function iteratively: simulate → run short chain → compute $\hat{p}_j$.

5. **Compute cppp and variance**
   - `compute_cppp(p_hat_obs, p_hat_cal)`
   - `transfer_ess_variance(delta_chain, p_hat_obs, p_hat_cal, m_tilde)`

---

<!--
Algorithm
**A. Run long chain**
1. Run long MCMC -> $\Delta$  chain
2. $\hat p_{\text{obs}} = M^{-1}\sum \mathbf{1}\{\Delta_i \le 0\}$

**B. Calibration (r short runs)**
For each replicate j:

  - Simulate $\tilde y_j$
  - Run short MCMC ($\tilde m_j$ iterations)
  - Compute $\hat p_j = \tilde m_j^{-1}\sum \mathbf{1}\{\tilde\Delta_{j,t}\le0\}$

**C. Estimate cppp**

$$
\widehat{\text{cppp}} = r^{-1}\sum_{j=1}^r
  \mathbf{1}\{\hat p_j \le \hat p_{\text{obs}}\}.
$$

**D. Variance via transfer ESS**
For each j:

  1. Let $q_j=\hat p_j$.
  2. Find threshold $q_j^{*}$ so empirical CDF of Δ at $q_j^*$= $q_j$.
  3. Form indicator $Z_i^{(j)}=\mathbf{1}\{\Delta_i\le q_j^*\}$.
  4. Estimate integrated autocorrelation $\tau_j$,
  5. Transfer ESS $\widetilde{\text{ESS}}_j=\tilde m_j/\tilde\tau_j$.
  6. $\widehat{\mathrm{Var}}(\hat{cppp})$


**E. Variance via two-piece decomposition (paper formulation)**

The cppp variance decomposes into two parts:

$$
\operatorname{Var}[\widehat{\text{cppp}}(y)]
=
\frac{1}{r}
\,\mathbb{E}_Y
\Big[
  \operatorname{Var}_{K|Y}(
    \mathbf{1}\{K \le \tilde m\,\widehat{\mathrm{ppp}}(y)\}\mid Y)
\Big]
+
\operatorname{Var}_Y
\Big[
  \mathbb{E}_{K|Y}(
    \mathbf{1}\{K \le \tilde m\,\widehat{\mathrm{ppp}}(y)\}\mid Y)
\Big].
$$
-->

##   Future extensions
- Support for other MCMC frameworks (`runCalibration_other()`).

---

## Example

```r
# define simulator for new datasets
new_data_fun <- make_data_sim_fun(...)

# create a discrepancy extractor
disc_fun <- make_col_disc_fun("ppp_column")

# or an offline discrepancy
disc_fun <- make_offline_disc_fun(control = list(
  new_data_fun = new_data_fun,
  discrepancy  = discrepancy))


# orchestrate calibration
res <- runCalibration(MCMC_samples_obs, MCMC_fun, new_data_fun, disc_fun,
                      num_reps = 50, control = list(seed = 123))

# compute cppp and variance
cppp_val <- compute_cppp(res$p_hat_obs, res$p_hat_cal)
var_out  <- transfer_ess_variance(res$delta_chain, res$p_hat_obs, res$p_hat_cal,
                                  m_tilde = res$m_tilde, c = 1.3)

summary(var_out)
```
