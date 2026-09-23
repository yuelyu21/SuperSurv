# Penalized Smooth Hazard Survival Learner

Fits an overall hazard model using
[`survPen::survPen()`](https://rdrr.io/pkg/survPen/man/survPen.html) and
returns its native survival probabilities. Only single-event,
right-censored outcomes and uniform observation weights are supported.

## Usage

``` r
surv.survPen(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  formula = NULL,
  baseline.df = 4L,
  n.legendre = 50L,
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

- formula:

  Optional one-sided hazard formula using feature names and
  `.supersurv_time`. For example,
  `~ smf(.supersurv_time, df = 4) + smf(x1, df = 4) + x2`, or
  `~ tensor(.supersurv_time, x1, df = c(4, 4)) + x2` for a time-varying
  effect. The `survPen` constructors `smf`, `tensor`, `tint`, and `rd`
  are available without attaching the backend package. Formula variables
  must come from the current training features or `.supersurv_time`.

- baseline.df:

  Baseline smooth degrees of freedom, at least 3. Used only when
  `formula` is NULL; the default adds linear terms for all features to
  `smf(.supersurv_time, df = baseline.df)`.

- n.legendre:

  Positive integer quadrature order for both fitting and survival
  prediction. Increase this to check integration accuracy.

- ...:

  Additional named fitting controls passed to
  [`survPen::survPen()`](https://rdrr.io/pkg/survPen/man/survPen.html),
  such as `lambda`, `method`, or `max.it.beta`.

## Value

A list with numeric survival matrix `pred` and fitted object `fit`.

## Details

Nonuniform observation weights are rejected because the backend does not
expose case-weighted fitting. The adapter does not implement net
survival, relative mortality, or left truncation. Survival curves that
increase beyond numerical tolerance cause an error; increase quadrature
accuracy rather than silently projecting a materially invalid curve.
Custom formulas must remain valid after feature screening; use
`screen.all` when explicitly naming features.

## Examples

``` r
if (requireNamespace("survPen", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  X <- dat[, "x1", drop = FALSE]
  fit <- surv.survPen(dat$duration, dat$event, X,
                      X[1:3, , drop = FALSE], c(50, 100), baseline.df = 3)
  dim(fit$pred)
}
#> [1] 3 2
```
