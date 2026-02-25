# Prediction function for GBM wrapper

Obtains predicted survivals from a fitted `surv.gbm` object. Uses a step
function to align the Breslow baseline hazard with the requested times.

## Usage

``` r
# S3 method for class 'surv.gbm'
predict(object, newX, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.gbm` object.

- newX:

  New covariate data.frame.

- new.times:

  Times at which to obtain predicted survivals.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newX` and columns correspond to the
evaluation times in `new.times`.
