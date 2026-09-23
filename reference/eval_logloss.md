# IPCW Log-Loss and Integrated IPCW Log-Loss

Evaluates the two-term inverse-probability-of-censoring weighted
log-loss. By default, censoring survival is estimated using marginal
reverse Kaplan-Meier. Conditional or externally estimated censoring
survival values can instead be supplied explicitly. Probability clipping
is applied only inside logarithms.

## Usage

``` r
eval_logloss(
  time,
  event,
  S_mat,
  times,
  tmin = min(times),
  tmax = max(times),
  ipcw_floor = 1e-06,
  eps = 1e-10,
  ipcw_cap = Inf,
  G_T_left = NULL,
  G_times = NULL
)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- S_mat:

  A numeric matrix of predicted survival probabilities (rows =
  observations, columns = time points).

- times:

  Numeric vector of evaluation times matching the columns of `S_mat`.

- tmin:

  Numeric. Lower bound for IBS integration. Defaults to `min(times)`.

- tmax:

  Numeric. Upper bound for IBS integration. Defaults to `max(times)`.

- ipcw_floor:

  Positive lower bound applied to censoring survival before inversion.

- eps:

  Probability clipping constant in `(0, 0.5)`.

- ipcw_cap:

  Largest IPCW contribution. Use `Inf` for no cap.

- G_T_left:

  Optional numeric vector containing subject-specific \\G(T_i-\mid
  X_i)\\ values. Supply together with `G_times`.

- G_times:

  Optional numeric matrix with the dimensions of `S_mat` containing
  \\G(t\mid X_i)\\, or a vector aligned with `times` for a common
  censoring survival curve. Supply together with `G_T_left`.

## Value

A list containing `logloss_scores`, `integrated_logloss`, `times`, and
`diagnostics`. Diagnostics report the effective denominator threshold,
intervention counts, weight quantiles, effective sample sizes, and the
fraction of weighted loss contributed by capped rows, separately for
failure and survivor terms.
