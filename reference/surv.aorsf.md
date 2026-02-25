# Wrapper for AORSF (Oblique Random Survival Forest)

Final Production Wrapper for AORSF (Tunable & Robust).

## Usage

``` r
surv.aorsf(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights,
  id,
  n_tree = 500,
  leaf_min_events = 5,
  mtry = NULL,
  ...
)
```

## Arguments

- time:

  Observed follow-up time; i.e. minimum of the event and censoring
  times.

- event:

  Observed event indicator; i.e, whether the follow-up time corresponds
  to an event or censoring.

- X:

  Training covariate data.frame.

- newX:

  Test covariate data.frame to use for prediction. Should have the same
  variable names and structure as `X`.

- new.times:

  Times at which to obtain the predicted survivals.

- obsWeights:

  Observation weights.

- id:

  Optional cluster/individual ID indicator.

- n_tree:

  Number of trees to grow (default: 500).

- leaf_min_events:

  Minimum number of events in a leaf node (default: 5).

- mtry:

  Number of predictors evaluated at each node.

- ...:

  Additional arguments passed to
  [`orsf`](https://docs.ropensci.org/aorsf/reference/orsf.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
