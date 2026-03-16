# Predict Method for Kaplan-Meier Wrapper

Obtains predicted survivals from a fitted `surv.km` object.

## Usage

``` r
# S3 method for class 'surv.km'
predict(object, newdata, new.times, ...)
```

## Arguments

- object:

  A fitted object of class `surv.km`.

- newdata:

  New covariate data.frame for which to obtain predictions (Ignored).

- new.times:

  Numeric vector of times at which to predict survival.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newdata` and columns correspond to
the evaluation times.

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.km(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

pred <- predict(fit$fit, newdata = newX, new.times = times)
dim(pred)
#> [1] 5 3
```
