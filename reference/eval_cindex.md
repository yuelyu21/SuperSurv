# Calculate Concordance Index (Harrell's or Uno's)

Calculate Concordance Index (Harrell's or Uno's)

## Usage

``` r
eval_cindex(time, event, S_mat, times, eval_time, method = "uno")
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- S_mat:

  A numeric matrix of predicted survival probabilities.

- times:

  Numeric vector of evaluation times matching the columns of `S_mat`.

- eval_time:

  Numeric. The specific time point at which to extract predictions.

- method:

  Character. Either "harrell" or "uno". Defaults to "uno".

## Value

A numeric value representing the chosen C-index.
