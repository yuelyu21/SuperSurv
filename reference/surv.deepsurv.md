# Experimental DeepSurv Learner

Fits a Cox partial-likelihood neural network through
[`survivalmodels::deepsurv()`](https://rdrr.io/pkg/survivalmodels/man/deepsurv.html).
This optional adapter requires a Python environment containing `torch`,
`torchtuples`, and `pycox` and currently supports only uniform
observation weights.

## Usage

``` r
surv.deepsurv(
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

- ...:

  Additional arguments passed to
  [`survivalmodels::deepsurv()`](https://rdrr.io/pkg/survivalmodels/man/deepsurv.html).

## Value

A list with numeric matrix `pred` and fitted object `fit`.

## Details

Fitted Python objects can be reused in the active R/Python session.
Plain [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) is not a
portable persistence format for these objects. Native survival curves
are mapped to requested times as right-continuous steps; see
[`surv.coxtime()`](https://yuelyu21.github.io/SuperSurv/reference/surv.coxtime.md)
for the boundary convention.

## Examples

``` r
if (interactive() && requireNamespace("survivalmodels", quietly = TRUE) &&
    requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_module_available("pycox")) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:60, ]
  X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
  fit <- surv.deepsurv(
    dat$duration, dat$event, X, X[1:4, , drop = FALSE],
    c(50, 100), epochs = 2, seed = 1
  )
  dim(fit$pred)
}
```
