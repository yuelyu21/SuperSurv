# Explain Predictions with Global SHAP (Kernel SHAP)

Explain Predictions with Global SHAP (Kernel SHAP)

## Usage

``` r
explain_kernel(
  model,
  X_explain,
  X_background,
  nsim = 20,
  only_best = FALSE,
  verbose = FALSE
)
```

## Arguments

- model:

  A fitted SuperSurv object OR a single wrapper output.

- X_explain:

  The dataset you want to explain (e.g., `X_test[1:10, ]`).

- X_background:

  The reference dataset for fastshap (e.g., `X_train[1:100, ]`).

- nsim:

  Number of simulations. Defaults to 20.

- only_best:

  Logical. If TRUE and model is SuperSurv, only explains the
  highest-weighted base learner.

- verbose:

  Logical; if `TRUE`, progress messages are shown.

## Value

A data.frame of class `c("explain", "data.frame")` containing the
calculated SHAP values. The columns correspond to the covariates in
`X_explain`.

## Examples

``` r
if (requireNamespace("fastshap", quietly = TRUE) &&
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
    X_explain = X[1:10, , drop = FALSE],
    X_background = X[11:40, , drop = FALSE],
    nsim = 5
  )

  dim(shap_values)
}
#> [1] 10  5
```
