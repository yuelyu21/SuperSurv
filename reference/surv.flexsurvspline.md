# Flexible Parametric Spline Survival Learner

Fits a weighted Royston-Parmar flexible parametric survival model using
[`flexsurv::flexsurvspline()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## Usage

``` r
surv.flexsurvspline(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  k = 1L,
  scale = "hazard",
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

- k:

  Number of internal spline knots.

- scale:

  Scale on which the flexible model is defined.

- ...:

  Additional arguments passed to
  [`flexsurv::flexsurvspline()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvspline.md).

## Value

A list with numeric matrix `pred` and fitted object `fit`.

## Examples

``` r
if (requireNamespace("flexsurv", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
  fit <- surv.flexsurvspline(
    dat$duration, dat$event, X, X[1:4, , drop = FALSE],
    c(50, 100), k = 1
  )
  dim(fit$pred)
}
#> [1] 4 2
```
