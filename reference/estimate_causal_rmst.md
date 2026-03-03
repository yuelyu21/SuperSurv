# Estimate Causal Restricted Mean Survival Time (RMST)

Calculates the causal treatment effect based on the difference in
Restricted Mean Survival Time (RMST) between treatment and control
groups up to a specific truncation time.

## Usage

``` r
estimate_causal_rmst(fit, data, trt_col, times, tau)
```

## Arguments

- fit:

  A fitted `SuperSurv` ensemble object.

- data:

  A `data.frame` containing the patient covariates and the treatment
  assignment.

- trt_col:

  Character string. The exact name of the binary treatment indicator
  column in `data` (e.g., "treatment").

- times:

  Numeric vector of time points matching the prediction grid.

- tau:

  Numeric. A single truncation time limit up to which the RMST will be
  calculated.

## Value

A numeric value representing the estimated causal RMST difference
(Treatment - Control).
