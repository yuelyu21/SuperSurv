# IPCW Brier Score and Integrated Brier Score (IBS)

Calculates the Inverse Probability of Censoring Weighted (IPCW) Brier
Score over a grid of times, and computes the Integrated Brier Score
(IBS) using trapezoidal integration.

## Usage

``` r
eval_brier(time, event, S_mat, times, tmin = min(times), tmax = max(times))
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

## Value

A list containing:

- `brier_scores`: A numeric vector of Brier scores at each time point.

- `ibs`: The Integrated Brier Score over the range
  `\code{tmin}, \code{tmax}`.

- `times`: The time grid used.
