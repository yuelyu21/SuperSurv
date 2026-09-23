# Flexible Parametric Distribution Survival Learner

Fits a weighted parametric survival model using
[`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
Native survival probabilities are used without risk-score calibration.

## Usage

``` r
surv.flexsurvreg(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  dist = "gengamma",
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

- dist:

  Distribution: `"gengamma"` (generalized gamma, default), `"gompertz"`,
  `"gamma"`, `"weibull"`, `"exp"`, `"lnorm"`, or `"llogis"`.

- ...:

  Additional named arguments to
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md),
  such as `inits`, `anc`, or optimizer `control`. Outcome, data, weight,
  truncation, and relative-survival arguments are managed or excluded by
  this adapter.

## Value

A list with numeric survival matrix `pred` and fitted object `fit`.

## Details

Only single-event, right-censored outcomes are supported. Covariates
enter the distribution's location parameter by default; `anc` can
specify covariates on ancillary parameters. Generalized gamma may
require suitable initial values, especially in small training folds.
Failed optimization is reported as an error rather than silently
accepted.

## Examples

``` r
if (requireNamespace("flexsurv", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  X <- dat[, "x1", drop = FALSE]
  fit <- surv.flexsurvreg(dat$duration, dat$event, X,
                         X[1:3, , drop = FALSE], c(50, 100), dist = "weibull")
  predict(fit$fit, X[1:3, , drop = FALSE], new.times = c(25, 50, 100))
}
#>           [,1]      [,2]      [,3]
#> [1,] 0.9331806 0.8539393 0.6973283
#> [2,] 0.9479166 0.8850395 0.7566717
#> [3,] 0.9227779 0.8323605 0.6577469
```
