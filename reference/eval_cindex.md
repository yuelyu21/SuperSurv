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

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:10, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.coxph(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

eval_cindex(
  time = dat$duration[1:10],
  event = dat$event[1:10],
  S_mat = fit$pred,
  times = times,
  eval_time = 100,
  method = "uno"
)
#> [1] 0.5898251
```
