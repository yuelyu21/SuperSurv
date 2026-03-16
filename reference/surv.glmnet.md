# Wrapper function for Penalized Cox Regression (GLMNET)

Final Production Wrapper for GLMNET (Tunable & Robust). Estimates a
penalized Cox model (Lasso, Ridge, or Elastic Net) with automatic lambda
selection. Uses the Breslow estimator with a step-function approach for
the baseline hazard.

## Usage

``` r
surv.glmnet(
  time,
  event,
  X,
  newdata,
  new.times,
  obsWeights,
  id,
  alpha = 1,
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

- newdata:

  Test covariate data.frame to use for prediction.

- new.times:

  Times at which to obtain the predicted survivals.

- obsWeights:

  Observation weights.

- id:

  Optional cluster/individual ID indicator.

- alpha:

  The elasticnet mixing parameter (0 = Ridge, 1 = Lasso). Default is 1.

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

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.glmnet(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    alpha = 1,
    nfolds = 3
  )

  dim(fit$pred)
}
#> [1] 5 3
```
