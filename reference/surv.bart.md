# Wrapper for BART (Bayesian Additive Regression Trees)

Final Production Wrapper for BART (Tunable & Robust). Uses the
[`mc.surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html) function.
Automatically reshapes the flat output vector into a survival matrix and
interpolates the predictions to the requested new.times.

## Usage

``` r
surv.bart(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  ntree = 10,
  ndpost = 30,
  nskip = 10,
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

  Times at which to obtain predicted survivals.

- obsWeights:

  Observation weights (Note: BART does not natively support weights).

- id:

  Optional cluster/individual ID indicator.

- ntree:

  Number of trees (default: 50).

- ndpost:

  Number of posterior draws (default: 1000).

- nskip:

  Number of burn-in draws (default: 250).

- ...:

  Additional arguments passed to
  [`mc.surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.

## Examples

``` r
if (.Platform$OS.type != "windows" &&
  requireNamespace("BART", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:20, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.bart(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    ntree = 3,
    ndpost = 5,
    nskip = 5
  )

  dim(fit$pred)
}
```
