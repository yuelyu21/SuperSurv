# Plot Predicted Survival Curves

Plot Predicted Survival Curves

## Usage

``` r
plot_predict(preds, eval_times, patient_idx = 1)
```

## Arguments

- preds:

  A list containing SuperSurv predictions OR a raw prediction matrix.

- eval_times:

  Numeric vector of times at which predictions were evaluated.

- patient_idx:

  Integer vector. Defaults to 1.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:120, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  eval_times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = eval_times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  preds <- predict(fit, newdata = X, new.times = eval_times)

  plot_predict(
    preds = preds,
    eval_times = eval_times,
    patient_idx = c(1, 2)
  )
}
```
