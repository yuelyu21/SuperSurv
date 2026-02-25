# Wrapper for Ridge Regression (Penalized Cox)

Final Production Wrapper for Ridge Regression (Tunable & Robust).
Estimates a penalized Cox model using a pure Ridge penalty (alpha = 0).

## Usage

``` r
surv.ridge(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights = NULL,
  id = NULL,
  nfolds = 10,
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

  Observation weights.

- id:

  Optional cluster/individual ID indicator.

- nfolds:

  Number of folds for internal cross-validation to select lambda.
  Default is 10.

- ...:

  Additional arguments passed to
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
