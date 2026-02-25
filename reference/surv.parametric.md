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
  newX,
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

- newX:

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
