# Super Learner for conditional survival functions

Orchestrates the cross-validation, metalearner optimization, and
prediction for an ensemble of survival base learners.

## Usage

``` r
SuperSurv(
  time,
  event,
  X,
  newX = NULL,
  new.times,
  event.SL.library,
  cens.SL.library,
  id = NULL,
  verbose = FALSE,
  control = list(),
  cvControl = list(),
  obsWeights = NULL,
  metalearner = "least_squares",
  selection = "ensemble",
  nFolds = 10,
  parallel = FALSE
)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data.frame.

- newX:

  Test covariate data.frame for prediction (defaults to X).

- new.times:

  Times at which to obtain predicted survivals.

- event.SL.library:

  Character vector of prediction algorithms for the event.

- cens.SL.library:

  Character vector of prediction algorithms for censoring.

- id:

  Cluster identification variable.

- verbose:

  Logical. If TRUE, prints progress messages.

- control:

  List of control parameters for the Super Learner.

- cvControl:

  List of control parameters for cross-validation.

- obsWeights:

  Observation weights.

- metalearner:

  Character string specifying the optimizer (e.g., "least_squares").

- selection:

  Character. Specifies how the meta-learner combines the base models.
  Use `"ensemble"` (default) to calculate a weighted average (convex
  combination) of the base learners. Use `"best"` to act as a Discrete
  Super Learner, which assigns a weight of 1.0 to the single model with
  the lowest cross-validated risk.

- nFolds:

  Number of cross-validation folds (default: 10).

- parallel:

  Logical. If TRUE, uses future.apply for parallel execution.

## Value

A list of class `SuperSurv` containing:

- `call`: The matched function call.

- `event.SL.predict`: Matrix of in-sample cross-validated survival
  predictions.

- `cens.SL.predict`: Matrix of in-sample cross-validated censoring
  predictions.

- `event.coef`: Numeric vector of optimized ensemble weights for the
  event.

- `cens.coef`: Numeric vector of optimized ensemble weights for
  censoring.

- `event.library.predict`: 3D array of cross-validated predictions from
  individual event learners.

- `event.libraryNames`: Data frame detailing the algorithms and
  screeners used.

- `event.fitLibrary`: List of the fitted base learner models (if
  `saveFitLibrary = TRUE`).

- `times`: The time grid used for evaluation.
