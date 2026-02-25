# Wrapper function for Random Survival Forests (RFSRC)

Final Production Wrapper for RFSRC (Tunable & Robust). Estimates a
survival random forest using
[`rfsrc`](https://www.randomforestsrc.org//reference/rfsrc.html).

## Usage

``` r
surv.rfsrc(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights = NULL,
  id = NULL,
  ntree = 1000,
  nodesize = 15,
  mtry = NULL,
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

  Times at which to obtain the predicted survivals.

- obsWeights:

  Observation weights.

- id:

  Currently ignored.

- ntree:

  Number of trees to grow (default: 1000).

- nodesize:

  Minimum number of deaths in terminal nodes (default: 15).

- mtry:

  Number of variables randomly selected as candidates for splitting a
  node.

- ...:

  Additional arguments passed to
  [`rfsrc`](https://www.randomforestsrc.org//reference/rfsrc.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
