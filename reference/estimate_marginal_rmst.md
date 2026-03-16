# Estimate Causal Restricted Mean Survival Time (RMST)

Calculates the causal treatment effect based on the difference in
Restricted Mean Survival Time (RMST) between treatment and control
groups up to a specific truncation time.

## Usage

``` r
estimate_marginal_rmst(fit, data, trt_col, times, tau, verbose = FALSE)
```

## Arguments

- fit:

  A fitted `SuperSurv` ensemble object.

- data:

  A `data.frame` containing the patient covariates and the treatment
  assignment.

- trt_col:

  Character string. The exact name of the binary treatment indicator
  column in `data` (e.g., "treatment").

- times:

  Numeric vector of time points matching the prediction grid.

- tau:

  Numeric. A single truncation time limit up to which the RMST will be
  calculated.

- verbose:

  Logical; if `TRUE`, progress messages are shown.

## Value

A numeric value representing the estimated causal RMST difference
(Treatment - Control).

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

rmst_res <- estimate_marginal_rmst(
  fit = fit,
  data = dat,
  trt_col = "x4",
  times = new.times,
  tau = 100
)
rmst_res$ATE_RMST
#> [1] -1.85712
```
