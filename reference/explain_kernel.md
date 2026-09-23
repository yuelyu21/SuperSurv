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
  verbose = FALSE,
  eval_time = NULL
)
```

## Arguments

- model:

  A fitted SuperSurv object OR a single wrapper output.

- X_explain:

  The dataset you want to explain (e.g., `X_test[1:10, ]`).

- X_background:

  Reference data defining the Kernel SHAP background distribution (for
  example, `X_train[1:100, ]`).

- nsim:

  Positive integer controlling the approximate coalition-sampling
  budget. It is converted to an even `m = 2 * nsim`; small feature sets
  are evaluated exactly by kernelshap. Defaults to 20.

- only_best:

  Logical. If TRUE and model is SuperSurv, only explains the
  highest-weighted base learner.

- verbose:

  Logical; if `TRUE`, progress messages are shown.

- eval_time:

  One finite, non-negative prediction time. Required explicitly so that
  all explanations target event probability `1 - S(eval_time)`.

## Value

A data.frame of class `c("explain", "data.frame")` containing the
calculated SHAP values. The columns correspond to the covariates in
`X_explain`. Attributes include `baseline`, `predictions`, `eval_time`,
`target = "event_probability"`, and backend convergence information.

## Details

The explained function uses the stored survival-prediction methods,
including screening and calibration. All positive ensemble weights are
used; small weights are not discarded. This replaces the earlier mixture
of native learner scores, so SHAP values from earlier releases are not
comparable. The background sample defines marginal, not conditional or
causal, SHAP values.

## Examples

``` r
if (FALSE) {
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
    nsim = 5, eval_time = 100
  )

  dim(shap_values)
}
```
