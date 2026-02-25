# Wrapper for standard Cox Proportional Hazards

Final Production Wrapper for CoxPH. Uses partial maximum likelihood and
the Breslow estimator.

## Usage

``` r
surv.coxph(time, event, X, newX, new.times, obsWeights, id, ...)
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

- ...:

  Additional arguments passed to
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
