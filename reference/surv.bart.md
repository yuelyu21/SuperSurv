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
  newX,
  new.times,
  obsWeights,
  id,
  ntree = 30,
  ndpost = 200,
  nskip = 100,
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
