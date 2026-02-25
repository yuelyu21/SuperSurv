# Predict Method for Kaplan-Meier Wrapper

Obtains predicted survivals from a fitted `surv.km` object.

## Usage

``` r
# S3 method for class 'surv.km'
predict(object, newX, new.times, ...)
```

## Arguments

- object:

  A fitted object of class `surv.km`.

- newX:

  New covariate data.frame for which to obtain predictions (Ignored).

- new.times:

  Numeric vector of times at which to predict survival.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newX` and columns correspond to the
evaluation times.
