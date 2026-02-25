# Explain Predictions with Global SHAP (Kernel SHAP)

Explain Predictions with Global SHAP (Kernel SHAP)

## Usage

``` r
explain_shap(model, X_explain, X_background, nsim = 20, only_best = FALSE)
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

## Value

A data.frame of class `c("explain", "data.frame")` containing the
calculated SHAP values. The columns correspond to the covariates in
`X_explain`.
