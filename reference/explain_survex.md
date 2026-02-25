# Create a Time-Dependent Survex Explainer

Bridges a fitted SuperSurv ensemble or a single base learner to the
`survex` package for Time-Dependent SHAP and Model Parts.

## Usage

``` r
explain_survex(model, data, y, times, label = NULL)
```

## Arguments

- model:

  A fitted SuperSurv object OR a single wrapper output.

- data:

  Covariate data for explanation (data.frame).

- y:

  The survival object (`Surv(time, event)`).

- times:

  The time grid for evaluation.

- label:

  Optional character string to name the explainer.

## Value

An explainer object of class `survex_explainer` created by
[`explain_survival`](https://modeloriented.github.io/survex/reference/explain_survival.html),
which can be passed to DALEX and survex functions for further model
diagnostics and plotting.
