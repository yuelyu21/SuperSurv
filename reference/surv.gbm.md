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
  newdata,
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

- newdata:

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

## Examples

``` r
if (requireNamespace("gbm", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.gbm(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    n.trees = 20,
    interaction.depth = 1,
    shrinkage = 0.05,
    cv.folds = 0,
    n.minobsinnode = 3
  )

  dim(fit$pred)
}
#> OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv_folds>1 when calling gbm usually results in improved predictive performance.
#> [1] 5 3
```
