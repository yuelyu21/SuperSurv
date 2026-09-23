# Evaluate survival predictions across models and times

Computes a common set of numerical benchmark results for a fitted
`SuperSurv` object or supported standalone learner. The censoring
distribution used by the Brier score and time-dependent AUC is estimated
marginally by Kaplan-Meier. For conditional censoring models,
resampling, inference, or formal model comparisons, use a specialist
evaluator such as
[`riskRegression::Score()`](https://rdrr.io/pkg/riskRegression/man/Score.html).

## Usage

``` r
eval_benchmark(
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

A list of class `"SuperSurv_benchmark"` containing `summary`, a
model-level table; `by_time`, a time-specific table; and the prediction
grid and risk horizon.
