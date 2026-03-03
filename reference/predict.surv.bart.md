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
