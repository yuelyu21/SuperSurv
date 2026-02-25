# Prediction function for Ranger wrapper

Obtains predicted survivals from a fitted `surv.ranger` object.

## Usage

``` r
# S3 method for class 'surv.ranger'
predict(object, newX, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.ranger` object.

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
