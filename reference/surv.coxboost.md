# Wrapper function for Component-Wise Boosting (CoxBoost)

Final Production Wrapper for CoxBoost (Tunable & Robust). Estimates a
Cox model via component-wise likelihood based boosting.

## Usage

``` r
surv.coxboost(
  time,
  event,
  X,
  newdata,
  new.times,
  obsWeights,
  id,
  stepno = 100,
  penalty = 100,
  ...
)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data.frame.

- newdata:

  Test covariate data.frame to use for prediction.

- new.times:

  Times at which to obtain the predicted survivals.

- obsWeights:

  Observation weights (Note: CoxBoost does not natively support weights,
  so these are ignored).

- id:

  Optional cluster/individual ID indicator.

- stepno:

  Number of boosting steps (default: 100).

- penalty:

  Penalty value for the update (default: 100).

- ...:

  Additional arguments passed to
  [`CoxBoost`](https://rdrr.io/pkg/CoxBoost/man/CoxBoost.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.

## Examples

``` r
if (requireNamespace("CoxBoost", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.coxboost(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    stepno = 10,
    penalty = 50
  )

  dim(fit$pred)
}
#> [1] 5 3
```
