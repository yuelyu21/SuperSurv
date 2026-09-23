# Control parameters for the SuperSurv Ensemble

Control parameters for the SuperSurv Ensemble

## Usage

``` r
SuperSurv.control(
  max.SL.iter = 20,
  event.t.grid = NULL,
  cens.t.grid = NULL,
  saveFitLibrary = TRUE,
  initWeightAlg = "surv.coxph",
  initWeight = "censoring",
  tol = 1e-05,
  traceIter = FALSE,
  ipcw.floor = 1e-04,
  ipcw.cap = 100,
  logloss.eps = 1e-10,
  optimizer.tol = 1e-08,
  optimizer.maxit = 10000L,
  truncation.warn.fraction = 0.05
)
```

## Arguments

- max.SL.iter:

  Maximum iterations for the iterative weighting algorithm. Default 20.

- event.t.grid:

  Optional time grid for event risk calculation.

- cens.t.grid:

  Optional time grid for censoring risk calculation.

- saveFitLibrary:

  Logical. If TRUE (default), saves models for future predictions.

- initWeightAlg:

  The learner used for the very first step of IPCW.

- initWeight:

  Whether to start by fitting "censoring" or "event" weights.

- tol:

  Positive convergence tolerance for the maximum change in the
  out-of-fold event and censoring ensemble predictions.

- traceIter:

  Logical. If TRUE, reports iteration diagnostics.

- ipcw.floor:

  Smallest censoring or event survival probability used in an IPCW
  denominator.

- ipcw.cap:

  Largest inverse-probability weight used in an IPCW loss.

- logloss.eps:

  Probability clipping constant used only while evaluating logarithms in
  the IPCW log-loss.

- optimizer.tol:

  Positive convergence tolerance for metalearner optimization.

- optimizer.maxit:

  Maximum number of metalearner optimization iterations.

- truncation.warn.fraction:

  Fraction of IPCW rows stabilized by flooring or truncation above which
  a warning is issued.

## Value

A list of control parameters.
