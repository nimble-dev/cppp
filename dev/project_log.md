# cpppPrototype — Development Log  

# 2025-11-17 
---

### Where to start next time
1. **Fix `new_data_fun` inside `runCalibrationNIMBLE()`:** I think there is some error in passing that 


### Where things stand
- `runCalibrationNIMBLE()` was modfied to:
  - Accept both uncompiled (`RmodelBaseClass`) and compiled (`CmodelBaseClass`) NIMBLE models.
  - Call the generic `runCalibration()` with `MCMC_fun`, `new_data_fun`, and `disc_fun`.

- Newcomb example:
  - Offline discrepancy `min_disc()` defined.
  - `disc_fun` built via `make_offline_disc_fun()` using a simple simulator `newcomb_newData()`.



# 2025-11-12 

### ✔ Unified Discrepancy Interface
- Standardized on: `disc_fun()` returns `list(obs = ..., sim = ...)`.
- Online mode with `make_col_disc_fun()`.
- Offline mode with `make_offline_disc_fun()`.
- Ensures one consistent discrepancy workflow across all backends.


### ✔ Refined `runCalibration()`
- Enforces correct `obs` / `sim` structure.
- Computes:
  - `PPP_obs`
  - `PPP_rep`
  - `CPPP`
- Added safety checks for lengths, structure, and row selection.

---

## **2025-11-12 — cPPP Computation Utilities**

### ✔ Implemented `computeCppp()`
Computes:
\[
\hat p_{cPPP} = \frac{1}{r} \sum I\{\text{pppCal}_j \le \text{pppObs}\}.
\]
Details:
- Validates scalar and vector inputs
- Handles empty calibration vectors gracefully

### ✔ Added tests for `computeCppp()`
Tests confirm:
- Correct numeric results for known inputs
- Monotonicity in `pppObs`
- Bounds in `[0,1]`
- Handling of empty calibration vectors

---

## **2025-11-12 — Roxygen, Math, and Documentation Fixes**

### ✔ Fixed Rd cross-reference warnings
- `[0,1]` incorrectly parsed as a link
- Replaced with `\eqn{[0,1]}` or escaped brackets `$begin:math:display$0,1$end:math:display$`

### ✔ Standardized mathematical expressions
Replaced `$...$` with:
- `\eqn{}` for inline math  
- `\deqn{}` for display math

---

## **Earlier — S3 Class and Base Infrastructure**

### ✔ Created S3 constructor: `new_cppresults()`
- Preserved original comments  
- Added `validate_cpppResult()`
- Checks:
  - `cppp ∈ [0,1]` when finite
  - `ppp` and `obs_ppp` numeric + within bounds  

### ✔ Constructor tests added
Ensures:
- Correct class structure  
- Proper lengths and numeric type  
- Bounds enforcement  

---

## **Earlier — Test Framework and Interface Scaffolding**

### ✔ testthat infrastructure set up
- Added `tests/testthat.R`
- Global helper for reproducible seeds

### ✔ Interface presence tests
Checked existence (not behavior) of:
- `computeCppp()`
- `transferAutocorrelation()`
- `runCalibration()`

These remain green even before full implementation.

---

## **Earlier — Package Setup & Cleanup**

### ✔ Cleaned and corrected `DESCRIPTION`
- Proper `Authors@R` formatting  
- Fixed license pointer  
- Removed duplicate `Suggests` entries  
- Added `URL`, `BugReports`, and `Config/testthat/edition: 3`  
- Added required newline at end of file  

### ✔ Added minimal `LICENSE`
Two-line file as required by R: