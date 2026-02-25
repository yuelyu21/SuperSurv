# Plot Survival Calibration Curve

Evaluates the calibration of the SuperSurv ensemble at a specific time
point by comparing predicted survival probabilities against observed
Kaplan-Meier estimates.

## Usage

``` r
plot_calibration(object, newdata, time, event, eval_time, bins = 5)
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

- eval_time:

  Numeric. A single time point at which to assess calibration.

- bins:

  Integer. The number of quantiles to group predictions into. Defaults
  to 5.

## Value

A ggplot object showing the calibration curve.
