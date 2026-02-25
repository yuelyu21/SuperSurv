# Prediction function for Cox regression wrapper

Obtains predicted survivals from a fitted `surv.coxph` object.

## Usage

``` r
# S3 method for class 'surv.coxph'
predict(object, newX, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.coxph` object.

- newX:

  New covariate data.frame for which to obtain predictions.

- new.times:

  Times at which to obtain predicted survivals.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newX` and columns correspond to the
evaluation times in `new.times`.
