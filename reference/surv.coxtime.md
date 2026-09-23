# Experimental Cox-Time Neural Survival Learner

Fits a neural Cox-Time model allowing time-dependent covariate effects
through
[`survivalmodels::coxtime()`](https://rdrr.io/pkg/survivalmodels/man/coxtime.html).
This optional adapter requires Python `torch`, `torchtuples`, and
`pycox` and supports only uniform observation weights.

## Usage

``` r
surv.coxtime(
  time,
  event,
  X,
  newdata = NULL,
  new.times,
  obsWeights = NULL,
  id = NULL,
  num_nodes = c(32L, 32L),
  activation = "relu",
  batch_norm = TRUE,
  dropout = NULL,
  epochs = 100L,
  batch_size = 128L,
  device = NULL,
  verbose = FALSE,
  seed = 1L,
  standardize_time = TRUE,
  ...
)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data frame.

- newdata:

  Covariate data frame used for prediction.

- new.times:

  Times at which survival probabilities are requested.

- obsWeights:

  Optional non-negative observation weights.

- id:

  Currently ignored.

- num_nodes:

  Positive integers giving hidden-layer sizes.

- activation:

  Neural-network activation name.

- batch_norm:

  Whether to use batch normalization.

- dropout:

  Optional dropout probability in `[0, 1)`.

- epochs:

  Number of training epochs.

- batch_size:

  Training and prediction batch size.

- device:

  Optional device passed to `survivalmodels`.

- verbose:

  Whether the Python backend should print training progress.

- seed:

  Positive integer used for Python, NumPy, and torch random-number
  generators.

- standardize_time:

  Whether to standardize outcome times using the training-only Cox-Time
  label transformation. Predictions are returned on the original time
  scale.

- ...:

  Additional named fitting arguments passed to
  [`survivalmodels::coxtime()`](https://rdrr.io/pkg/survivalmodels/man/coxtime.html),
  such as `frac` for an internal training-fold validation split or
  `early_stopping`. Outcomes and features are managed by the adapter.

## Value

A list with numeric survival matrix `pred` and fitted object `fit`.

## Details

Native survival curves are evaluated as right-continuous steps on the
backend's time grid, using survival one before the first grid point and
the final available value beyond the last point. This is not a claim of
reliable extrapolation beyond observed follow-up. Python-backed fitted
objects are intended for reuse within the active R/Python session; plain
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) is not a portable
persistence format for Python objects.

## Examples

``` r
if (interactive() && requireNamespace("survivalmodels", quietly = TRUE) &&
    requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_module_available("pycox")) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:60, ]
  X <- dat[, "x1", drop = FALSE]
  fit <- surv.coxtime(dat$duration, dat$event, X,
                      X[1:3, , drop = FALSE], c(50, 100),
                      epochs = 2, batch_norm = FALSE, device = "cpu")
  dim(fit$pred)
}
```
