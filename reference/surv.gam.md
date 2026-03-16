# Wrapper for Generalized Additive Cox Regression (GAM)

Final Production Wrapper for GAM (Tunable & Robust). Uses
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) to fit an additive
combination of smooth and linear functions.

## Usage

``` r
surv.gam(time, event, X, newdata, new.times, obsWeights, id, cts.num = 5, ...)
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

  Times at which to obtain predicted survivals.

- obsWeights:

  Observation weights (Note: Ignored, as mgcv uses weights for the event
  indicator).

- id:

  Optional cluster/individual ID indicator.

- cts.num:

  Cutoff of unique values at which a numeric covariate receives a smooth
  term (s).

- ...:

  Additional arguments passed to
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.

## Examples

``` r
if (requireNamespace("mgcv", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.gam(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    cts.num = 5
  )

  dim(fit$pred)
}
#> [1] 5 3
```
