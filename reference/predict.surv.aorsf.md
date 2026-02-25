# Prediction function for AORSF

Obtains predicted survivals from a fitted `surv.aorsf` object. Uses the
native `aorsf` prediction engine to calculate survival directly at
requested times.

## Usage

``` r
# S3 method for class 'surv.aorsf'
predict(object, newX, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.aorsf` object.

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
