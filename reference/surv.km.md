# Kaplan-Meier Prediction Algorithm

This prediction algorithm ignores all covariates and computes the
marginal Kaplan-Meier survival estimator using the
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) function.

## Usage

``` r
surv.km(time, event, X, newX, new.times, obsWeights, id, ...)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- X:

  Training covariate data.frame (Ignored by KM).

- newX:

  Test covariate data.frame to use for prediction.

- new.times:

  Numeric vector of times at which to predict survival.

- obsWeights:

  Numeric vector of observation weights.

- id:

  Optional vector indicating subject/cluster identities.

- ...:

  Additional ignored arguments.

## Value

A list containing:

- `fit`: A list containing the fitted
  [`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) object.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at `new.times`.
