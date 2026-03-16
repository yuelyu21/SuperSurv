# Universal Parametric Survival Wrapper

Final Production Wrapper for AFT Models (Weibull, Exponential,
LogNormal, LogLogistic). Replaces individual wrappers with one robust,
vectorized function.

## Usage

``` r
surv.parametric(
  time,
  event,
  X,
  newdata,
  new.times,
  obsWeights,
  id,
  dist = "weibull",
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

  Observation weights.

- id:

  Optional cluster/individual ID indicator.

- dist:

  Distribution for the AFT model (default: "weibull").

- ...:

  Additional arguments passed to
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.parametric(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL,
  dist = "weibull"
)

dim(fit$pred)
#> [1] 5 3
```
