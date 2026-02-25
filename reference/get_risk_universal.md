# Universal Risk Wrapper for SHAP

Ensures "High Value = High Risk" orientation for ALL base learners and
the SuperSurv ensemble, allowing unified SHAP value calculations.

## Usage

``` r
get_risk_universal(object, newdata)
```

## Arguments

- object:

  A fitted model object (e.g., from a single wrapper or a SuperSurv
  ensemble).

- newdata:

  A data.frame of new covariates to predict on.

## Value

A numeric vector of risk scores of the same length as the number of rows
in `newdata`, where higher values consistently indicate a higher risk of
the event.
