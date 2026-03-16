# Prediction function for BART

Obtains predicted survivals from a fitted `surv.bart` object. Manually
expands the time grid to bypass the 10-column expectation.

## Usage

``` r
# S3 method for class 'surv.bart'
predict(object, newdata, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.bart` object.

- newdata:

  New covariate data.frame to predict on.

- new.times:

  Times to predict.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newdata` and columns correspond to
the evaluation times in `new.times`.

## Examples

``` r
if (.Platform$OS.type != "windows" &&
  requireNamespace("BART", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:20, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.bart(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    ntree = 3,
    ndpost = 5,
    nskip = 5
  )

  pred <- fit$pred
  dim(pred)
}
```
