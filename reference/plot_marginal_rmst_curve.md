# Plot Adjusted Marginal RMST Contrast Over Time

Generates a curve showing how the adjusted marginal restricted mean
survival time (RMST) contrast evolves across a sequence of restriction
times.

## Usage

``` r
plot_marginal_rmst_curve(
  fit,
  data,
  trt_col,
  times,
  tau_seq,
  inference = FALSE,
  B = 200,
  seed = NULL,
  ci_level = 0.95
)
```

## Arguments

- fit:

  A fitted `SuperSurv` ensemble object.

- data:

  A `data.frame` containing the covariates and the binary grouping
  variable.

- trt_col:

  Character string. The exact name of the binary grouping variable in
  `data`.

- times:

  Numeric vector of time points matching the prediction grid.

- tau_seq:

  Numeric vector. A sequence of restriction times (`tau`) to evaluate
  and plot.

- inference:

  Deprecated compatibility argument. Must remain `FALSE`.

- B, seed, ci_level:

  Deprecated compatibility arguments; ignored.

## Value

A `ggplot` object visualizing the adjusted marginal RMST contrast curve.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
data("metabric", package = "SuperSurv")
dat <- metabric[1:80, ]
x_cols <- grep("^x", names(dat), value = TRUE)[1:5]
X <- dat[, x_cols, drop = FALSE]
new.times <- seq(20, 120, by = 20)

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X,
  new.times = new.times,
  event.library = c("surv.coxph"),
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
}
#> Adjusted Delta RMST at tau = 40: -0.158 time units
#> Adjusted Delta RMST at tau = 60: -0.694 time units
#> Adjusted Delta RMST at tau = 80: -1.643 time units
#> Adjusted Delta RMST at tau = 100: -2.728 time units
#> Adjusted Delta RMST at tau = 120: -4.056 time units

```
