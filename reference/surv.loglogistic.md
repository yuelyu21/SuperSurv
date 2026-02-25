# Parametric Survival Prediction Wrapper (Log-Logistic)

Parametric Survival Prediction Wrapper (Log-Logistic)

## Usage

``` r
surv.loglogistic(time, event, X, newX, new.times, obsWeights, id, ...)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data.frame.

- newX:

  Test covariate data.frame to use for prediction.

- new.times:

  Times at which to obtain the predicted survivals.

- obsWeights:

  Observation weights.

- id:

  Cluster identification variable.

- ...:

  Additional ignored arguments.

## Value

A list containing the fitted model and predictions.
