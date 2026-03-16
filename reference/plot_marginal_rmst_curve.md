# Plot Causal RMST Difference Over Time

Generates a curve showing how the causal Restricted Mean Survival Time
(RMST) difference between treatment groups evolves across a sequence of
different truncation times.

## Usage

``` r
plot_marginal_rmst_curve(fit, data, trt_col, times, tau_seq)
```

## Arguments

- fit:

  A fitted `SuperSurv` ensemble object.

- data:

  A `data.frame` containing the patient covariates and the treatment
  assignment.

- trt_col:

  Character string. The exact name of the binary treatment indicator
  column in `data`.

- times:

  Numeric vector of time points matching the prediction grid.

- tau_seq:

  Numeric vector. A sequence of truncation times (`tau`) to evaluate and
  plot.

## Value

A `ggplot` object visualizing the causal RMST difference curve.

## Examples

``` r
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
  event.library = c("surv.coxph", "surv.glmnet"),
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE)
)

tau_grid <- seq(40, 120, by = 20)
plot_marginal_rmst_curve(
  fit = fit,
  data = dat,
  trt_col = "x4",
  times = new.times,
  tau_seq = tau_grid
)
```
