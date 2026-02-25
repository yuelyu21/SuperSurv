# Survival Probability Heatmap

Survival Probability Heatmap

## Usage

``` r
plot_survival_heatmap(object, newdata, times)
```

## Arguments

- object:

  A fitted SuperSurv object

- newdata:

  Test covariates (e.g., `X_te[1:50, ]`)

- times:

  The time grid to visualize

## Value

A `ggplot` object visualizing the SHAP values.
