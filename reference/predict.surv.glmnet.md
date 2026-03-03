# Prediction function for GLMNET wrapper

Obtains predicted survivals from a fitted `surv.glmnet` object.

## Usage

``` r
# S3 method for class 'surv.glmnet'
predict(object, newdata, new.times, ...)
```

## Arguments

- object:

  Fitted `surv.glmnet` object.

- newdata:

  New covariate data.frame.

- new.times:

  Times at which to obtain predicted survivals.

- ...:

  Additional ignored arguments.

## Value

A numeric matrix of predicted survival probabilities, where rows
correspond to the observations in `newdata` and columns correspond to
the evaluation times in `new.times`.
