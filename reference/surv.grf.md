# Generalized Random Survival Forest Learner

Fits an honest generalized random survival forest using
[`grf::survival_forest()`](https://rdrr.io/pkg/grf/man/survival_forest.html)
and returns conditional survival curves.

## Usage

``` r
surv.grf(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  num.trees = 1000L,
  mtry = NULL,
  min.node.size = 15L,
  honesty = TRUE,
  prediction.type = c("Kaplan-Meier", "Nelson-Aalen"),
  seed = 1L,
  ...
)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data frame.

- newdata:

  Covariate data frame used for prediction.

- new.times:

  Times at which survival probabilities are requested.

- obsWeights:

  Optional non-negative observation weights.

- id:

  Currently ignored.

- num.trees:

  Number of trees.

- mtry:

  Number of candidate variables considered at each split.

- min.node.size:

  Minimum terminal-node size.

- honesty:

  Whether to use honest sample splitting.

- prediction.type:

  Either `"Kaplan-Meier"` or `"Nelson-Aalen"`.

- seed:

  Integer random seed passed to `grf`.

- ...:

  Additional arguments passed to
  [`grf::survival_forest()`](https://rdrr.io/pkg/grf/man/survival_forest.html).

## Value

A list with numeric matrix `pred` and fitted object `fit`.

## Examples

``` r
if (requireNamespace("grf", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
  fit <- surv.grf(
    dat$duration, dat$event, X, X[1:4, , drop = FALSE],
    c(50, 100), num.trees = 50, seed = 1
  )
  dim(fit$pred)
}
#> [1] 4 2
```
