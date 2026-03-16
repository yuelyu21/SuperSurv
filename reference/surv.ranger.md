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
  newdata,
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

- newdata:

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

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.ranger(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    num.trees = 10,
    min.node.size = 3
  )

  dim(fit$pred)
}
#> [1] 5 3
```
