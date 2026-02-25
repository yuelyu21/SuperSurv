# Wrapper for Survival Regression Trees (rpart)

Final Production Wrapper for single decision trees. Uses
[`rpart`](https://rdrr.io/pkg/rpart/man/rpart.html) with method="exp"
and calculates survival probabilities using the Breslow estimator.

## Usage

``` r
surv.rpart(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights,
  id,
  cp = 0.01,
  minsplit = 20,
  maxdepth = 30,
  ...
)
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

  Times at which to obtain predicted survivals.

- obsWeights:

  Observation weights.

- id:

  Optional cluster/individual ID indicator.

- cp:

  Complexity parameter (default: 0.01).

- minsplit:

  Minimum number of observations to attempt a split (default: 20).

- maxdepth:

  Maximum depth of any node of the final tree (default: 30).

- ...:

  Additional arguments passed to
  [`rpart.control`](https://rdrr.io/pkg/rpart/man/rpart.control.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
