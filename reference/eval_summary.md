# Evaluate SuperSurv predictions on test data

Computes the integrated Brier score (IBS), Uno C-index, and integrated
area under the curve (iAUC) for the SuperSurv ensemble and all
individual base learners.

## Usage

``` r
eval_summary(
  object,
  newdata,
  time,
  event,
  eval_times,
  risk_time = stats::median(eval_times),
  verbose = FALSE
)
```

## Arguments

- object:

  A fitted `SuperSurv` object.

- newdata:

  A data.frame of test covariates.

- time:

  Numeric vector of observed follow-up times for the test set.

- event:

  Numeric vector of event indicators for the test set.

- eval_times:

  Numeric vector of times at which to evaluate survival predictions.

- risk_time:

  Numeric. The specific time horizon used when extracting risk scores
  for Uno C-index. Defaults to the median of `eval_times`.

- verbose:

  Logical; if `TRUE`, progress messages are shown.

## Value

An object of class `"SuperSurv_eval"` containing benchmark metrics for
the ensemble and base learners.

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:80, ]
x_cols <- grep("^x", names(dat))[1:3]

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = dat[, x_cols, drop = FALSE],
  new.times = seq(50, 200, by = 50),
  event.library = c("surv.coxph", "surv.km"),
  cens.library = c("surv.coxph", "surv.km")
)

res <- eval_summary(
  object = fit,
  newdata = dat[, x_cols, drop = FALSE],
  time = dat$duration,
  event = dat$event,
  eval_times = seq(50, 200, by = 50)
)

res
#>                   Model    IBS  Uno_C   iAUC
#> 1    SuperSurv_Ensemble 0.2272 0.5713 0.5988
#> 2 surv.coxph_screen.all 0.2172 0.5713 0.5988
#> 3    surv.km_screen.all 0.2306 0.5000 0.5000
```
