# Plot Predicted RMST vs. Observed Survival Times

Evaluates the calibration of the causal RMST estimator by plotting the
model's predicted RMST for each patient against their actual observed
follow-up time.

## Usage

``` r
plot_rmst_vs_obs(fit, data, time_col, event_col, times, tau)
```

## Arguments

- fit:

  A fitted `SuperSurv` ensemble object.

- data:

  A `data.frame` containing the patient covariates, times, and events.

- time_col:

  Character string. The exact name of the observed follow-up time column
  in `data`.

- event_col:

  Character string. The exact name of the event indicator column in
  `data` (e.g., 1 for event, 0 for censored).

- times:

  Numeric vector of time points matching the prediction grid.

- tau:

  Numeric. A single truncation time limit up to which the RMST is
  calculated.

## Value

A `ggplot` object comparing predicted RMST to observed outcomes.
