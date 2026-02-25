# Calculate Baseline Survival using Breslow Estimator

Calculate Baseline Survival using Breslow Estimator

## Usage

``` r
safe_breslow_step(time, event, risk_score, new.times)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (0=censored, 1=event).

- risk_score:

  Numeric vector of risk scores.

- new.times:

  Numeric vector of times at which to evaluate the baseline hazard.

## Value

A numeric vector representing the baseline survival probabilities
evaluated at the specified `new.times`.
