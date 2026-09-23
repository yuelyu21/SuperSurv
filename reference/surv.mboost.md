# Component-Wise Cox Boosting Learner

Fits a weighted component-wise Cox proportional-hazards model using
[`mboost::glmboost()`](https://rdrr.io/pkg/mboost/man/glmboost.html) and
converts its risk score to survival probabilities using SuperSurv's
weighted baseline-hazard calibration.

## Usage

``` r
surv.mboost(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  mstop = 100L,
  nu = 0.1,
  center = FALSE,
  ties = c("breslow", "efron"),
  survival_transform = c("exponential", "product_limit"),
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

- mstop:

  Number of boosting iterations.

- nu:

  Boosting step size.

- center:

  Whether to center component-wise base learners.

- ties:

  Tied-event approximation for risk-score calibration.

- survival_transform:

  Transformation from calibrated hazard increments to survival
  probabilities.

- ...:

  Additional arguments passed to
  [`mboost::glmboost()`](https://rdrr.io/pkg/mboost/man/glmboost.html).

## Value

A list with numeric matrix `pred` and fitted object `fit`.

## Examples

``` r
if (requireNamespace("mboost", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
  fit <- surv.mboost(
    dat$duration, dat$event, X, X[1:4, , drop = FALSE],
    c(50, 100), mstop = 20
  )
  dim(fit$pred)
}
#> [1] 4 2
```
