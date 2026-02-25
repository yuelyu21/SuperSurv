# Control parameters for the SuperSurv Ensemble

Control parameters for the SuperSurv Ensemble

## Usage

``` r
SuperSurv.control(
  max.SL.iter = 20,
  event.t.grid = NULL,
  cens.t.grid = NULL,
  saveFitLibrary = TRUE,
  initWeightAlg = "surv.coxph",
  initWeight = "censoring",
  ...
)
```

## Arguments

- max.SL.iter:

  Maximum iterations for the iterative weighting algorithm. Default 20.

- event.t.grid:

  Optional time grid for event risk calculation.

- cens.t.grid:

  Optional time grid for censoring risk calculation.

- saveFitLibrary:

  Logical. If TRUE (default), saves models for future predictions.

- initWeightAlg:

  The learner used for the very first step of IPCW.

- initWeight:

  Whether to start by fitting "censoring" or "event" weights.

- ...:

  Additional ignored arguments.

## Value

A list of control parameters.
