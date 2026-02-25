# Plot Global Feature Importance for SuperSurv

Plot Global Feature Importance for SuperSurv

## Usage

``` r
plot_global_importance(
  shap_values,
  title = "SuperSurv: Ensemble Feature Importance",
  top_n = 10
)
```

## Arguments

- shap_values:

  The output from explain_shap.SuperSurv

- title:

  Plot title.

- top_n:

  Number of features to show (default 10)

## Value

A `ggplot` object visualizing the SHAP values.
