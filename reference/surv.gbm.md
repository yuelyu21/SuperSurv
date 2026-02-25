# Wrapper function for Gradient Boosting (GBM) prediction algorithm

Final Production Wrapper for GBM (Tunable & Robust). Estimates a Cox
proportional hazards model via gradient boosting. Uses the Breslow
estimator with a step-function approach for the baseline hazard.
Includes internal safeguards against C++ crashes and small
cross-validation folds.

## Usage

``` r
surv.gbm(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights,
  id,
  n.trees = 1000,
  interaction.depth = 2,
  shrinkage = 0.01,
  cv.folds = 5,
  n.minobsinnode = 10,
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

- n.trees:

  Integer specifying the total number of trees to fit (default: 1000).

- interaction.depth:

  Maximum depth of variable interactions (default: 2).

- shrinkage:

  A shrinkage parameter applied to each tree (default: 0.01).

- cv.folds:

  Number of cross-validation folds to perform internally for optimal
  tree selection (default: 5).

- n.minobsinnode:

  Minimum number of observations in the trees terminal nodes (default:
  10).

- ...:

  Additional arguments passed to
  [`gbm`](https://rdrr.io/pkg/gbm/man/gbm.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
