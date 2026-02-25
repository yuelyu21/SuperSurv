# Wrapper function for Component-Wise Boosting (CoxBoost)

Final Production Wrapper for CoxBoost (Tunable & Robust). Estimates a
Cox model via component-wise likelihood based boosting.

## Usage

``` r
surv.coxboost(
  time,
  event,
  X,
  newX,
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

- newX:

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
