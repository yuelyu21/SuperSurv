# Beeswarm Summary Plot for SuperSurv SHAP

Beeswarm Summary Plot for SuperSurv SHAP

## Usage

``` r
plot_beeswarm(shap_values, data, top_n = 10)
```

## Arguments

- shap_values:

  The output from
  [`explain_kernel()`](https://yuelyu21.github.io/SuperSurv/reference/explain_kernel.md).

- data:

  The covariate data used (X_explain)

- top_n:

  Number of features to display

## Value

A `ggplot` object visualizing the SHAP values.

## Examples

``` r
if (requireNamespace("fastshap", quietly = TRUE) &&
    requireNamespace("ggforce", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  shap_values <- explain_kernel(
    model = fit,
    X_explain = X[1:20, , drop = FALSE],
    X_background = X[21:50, , drop = FALSE],
    nsim = 5
  )

  plot_beeswarm(
    shap_values = shap_values,
    data = X[1:20, , drop = FALSE],
    top_n = 5
  )
}
```
