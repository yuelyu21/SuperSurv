# Wrapper function for Ranger Random Survival Forest

Final Production Wrapper for Ranger (Tunable & Fast). Uses the
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) C++
implementation to estimate survival curves.

## Usage

``` r
surv.ranger(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights,
  id,
  num.trees = 500,
  mtry = NULL,
  min.node.size = NULL,
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

- num.trees:

  Number of trees (default: 500).

- mtry:

  Number of variables to split at each node. Defaults to `sqrt(p)`.

- min.node.size:

  Minimum node size (default: 15 for survival).

- ...:

  Additional arguments passed to
  [`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
