# Predict method for SuperSurv fits

Obtains predicted survival probabilities from a fitted SuperSurv
ensemble.

## Usage

``` r
# S3 method for class 'SuperSurv'
predict(object, newdata, new.times, onlySL = FALSE, threshold = 1e-04, ...)
```

## Arguments

- object:

  A fitted object of class `SuperSurv`.

- newdata:

  A data.frame of new covariate values.

- new.times:

  A numeric vector of times at which to predict survival.

- onlySL:

  Logical. If TRUE, only uses models with weights \> threshold.

- threshold:

  Numeric. The weight threshold for onlySL.

- ...:

  Additional ignored arguments.

## Value

A list containing:

- `event.predict`: A numeric matrix (rows = observations, columns =
  times) of the final predicted survival probabilities from the
  ensemble.

- `event.library.predict`: A 3D numeric array (observations x times x
  models) containing the individual survival predictions from each base
  learner.

- `cens.predict`: A numeric matrix of the predicted censoring
  probabilities.

- `cens.library.predict`: A 3D numeric array of the individual censoring
  predictions.

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:10, , drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  preds <- predict(
    object = fit,
    newdata = newX,
    new.times = new.times
  )

  dim(preds$event.predict)
}
#> [1] 10  6
```
