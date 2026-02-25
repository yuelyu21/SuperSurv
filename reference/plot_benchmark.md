# Plot Longitudinal Benchmark Metrics

Generates time-dependent performance curves comparing the SuperSurv
ensemble against its base learners.

## Usage

``` r
plot_benchmark(
  object,
  newdata,
  time,
  event,
  eval_times,
  metrics = c("brier", "auc", "cindex")
)
```

## Arguments

- object:

  A fitted SuperSurv object.

- newdata:

  A data.frame of test covariates.

- time:

  Numeric vector of observed follow-up times for the test set.

- event:

  Numeric vector of event indicators for the test set.

- eval_times:

  Numeric vector of times at which to evaluate predictions.

- metrics:

  Character vector specifying which plots to return. Options: "brier",
  "auc", "cindex". Defaults to all three.

## Value

A combined patchwork ggplot object, or a single ggplot if only one
metric is selected.
