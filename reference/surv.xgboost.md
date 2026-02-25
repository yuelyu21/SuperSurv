# Wrapper for XGBoost (Robust CV-Tuned + Safe Prediction)

Estimates a Cox proportional hazards model via XGBoost. Incorporates
safe Breslow hazard calculation and matrix alignment to prevent C++
crashes.

## Usage

``` r
surv.xgboost(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights,
  id,
  nrounds = 1000,
  early_stopping_rounds = 10,
  eta = 0.05,
  max_depth = 2,
  min_child_weight = 5,
  lambda = 10,
  subsample = 0.7,
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

- nrounds:

  Max number of boosting iterations (default: 1000).

- early_stopping_rounds:

  Rounds with no improvement to trigger early stopping (default: 10).

- eta:

  Learning rate (default: 0.05).

- max_depth:

  Maximum tree depth (default: 2).

- min_child_weight:

  Minimum sum of instance weight in a child (default: 5).

- lambda:

  L2 regularization term on weights (default: 10).

- subsample:

  Subsample ratio of the training instances (default: 0.7).

- ...:

  Additional arguments passed to
  [`xgb.train`](https://rdrr.io/pkg/xgboost/man/xgb.train.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
