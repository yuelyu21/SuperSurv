# Prediction function for GBM wrapper

Obtains predicted survivals from a fitted `surv.gbm` object. Uses a step
function to align the Breslow baseline hazard with the requested times.

## Usage

``` r
# S3 method for class 'surv.gbm'
predict(object, newdata, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.gbm` object.

- newdata:

  New covariate data.frame.

- new.times:

  Times at which to obtain predicted survivals.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newdata` and columns correspond to
the evaluation times in `new.times`.

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

  pred <- predict(fit$fit, newdata = newX, new.times = times)
  dim(pred)
}
#> OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv_folds>1 when calling gbm usually results in improved predictive performance.
#> [1] 5 3
```
