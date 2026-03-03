# Plot Causal RMST Difference Over Time

Generates a curve showing how the causal Restricted Mean Survival Time
(RMST) difference between treatment groups evolves across a sequence of
different truncation times.

## Usage

``` r
plot_causal_rmst_curve(fit, data, trt_col, times, tau_seq)
```

## Arguments

- fit:

  A fitted `SuperSurv` ensemble object.

- data:

  A `data.frame` containing the patient covariates and the treatment
  assignment.

- trt_col:

  Character string. The exact name of the binary treatment indicator
  column in `data`.

- times:

  Numeric vector of time points matching the prediction grid.

- tau_seq:

  Numeric vector. A sequence of truncation times (`tau`) to evaluate and
  plot.

## Value

A `ggplot` object visualizing the causal RMST difference curve.
