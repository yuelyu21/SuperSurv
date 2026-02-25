# Prediction function for Universal Parametric Wrapper

Obtains predicted survivals from a fitted `surv.parametric` object using
exact closed-form equations.

## Usage

``` r
# S3 method for class 'surv.parametric'
predict(object, newX, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.parametric` object.

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
