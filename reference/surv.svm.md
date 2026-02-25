# Wrapper for Survival Support Vector Machine (survivalsvm)

Final Production Wrapper for SVM (Tunable & Robust). Estimates a
survival SVM and calibrates the raw utility scores into survival
probabilities using a univariate Cox proportional hazards model.

## Usage

``` r
surv.svm(
  time,
  event,
  X,
  newX,
  new.times,
  obsWeights,
  id,
  gamma.mu = 0.1,
  type = "vanbelle2",
  kernel = "lin_kernel",
  opt.meth = "quadprog",
  ...
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

  Test covariate data.frame to use for prediction.

- new.times:

  Times at which to obtain the predicted survivals.

- obsWeights:

  Observation weights.

- id:

  Optional cluster/individual ID indicator.

- gamma.mu:

  Regularization parameter for the SVM (default: 0.1).

- type:

  Type of SVM implementation (default: "vanbelle2").

- kernel:

  Kernel type for the SVM (default: "lin_kernel").

- opt.meth:

  Optimization method (default: "quadprog").

- ...:

  Additional arguments passed to
  [`survivalsvm`](https://rdrr.io/pkg/survivalsvm/man/survivalsvm.html).

## Value

A list containing:

- `fit`: The fitted model object (e.g., the raw `coxph` or `xgb.Booster`
  object). If the model fails to fit, this may be an object of class
  `try-error`.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at the specified `new.times` grid.
