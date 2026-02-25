# Evaluate SuperSurv Predictions on Test Data

Computes the Integrated Brier Score (IBS), Uno's C-index, and Integrated
AUC (iAUC) for the SuperSurv ensemble and all individual base learners.

## Usage

``` r
eval_summary(
  object,
  newdata,
  time,
  event,
  eval_times,
  risk_time = median(eval_times)
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

  Numeric. The specific time horizon to use when extracting risk scores
  for Uno's C-index. Defaults to the median of `eval_times`.

## Value

A data.frame containing the benchmark metrics for all models.
