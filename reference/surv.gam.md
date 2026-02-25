# Wrapper for Generalized Additive Cox Regression (GAM)

Final Production Wrapper for GAM (Tunable & Robust). Uses
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) to fit an additive
combination of smooth and linear functions.

## Usage

``` r
surv.gam(time, event, X, newX, new.times, obsWeights, id, cts.num = 5, ...)
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
