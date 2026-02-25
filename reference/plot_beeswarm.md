# Beeswarm Summary Plot for SuperSurv SHAP

Beeswarm Summary Plot for SuperSurv SHAP

## Usage

``` r
plot_beeswarm(shap_values, data, top_n = 10)
```

## Arguments

- shap_values:

  The output from explain_shap.SuperSurv

- data:

  The covariate data used (X_explain)

- top_n:

  Number of features to display

## Value

A `ggplot` object visualizing the SHAP values.
