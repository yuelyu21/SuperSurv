# Plot SHAP Dependence for SuperSurv

Plot SHAP Dependence for SuperSurv

## Usage

``` r
plot_dependence(shap_values, data, feature_name, title = NULL)
```

## Arguments

- shap_values:

  The output from explain_shap.SuperSurv

- data:

  The original covariate data used for the explanation (X_explain)

- feature_name:

  String name of the column to plot

- title:

  Optional custom title.

## Value

A `ggplot` object visualizing the SHAP values.
