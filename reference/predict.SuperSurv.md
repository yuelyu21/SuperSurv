# Predict method for SuperSurv fits

Obtains predicted survival probabilities from a fitted SuperSurv
ensemble.

## Usage

``` r
# S3 method for class 'SuperSurv'
predict(
  object,
  newdata,
  new.times,
  type = c("both", "event", "censoring"),
  onlySL = FALSE,
  threshold = 1e-04,
  ...
)
```

## Arguments

- object:

  A fitted object of class `SuperSurv`.

- newdata:

  A data.frame of new covariate values.

- new.times:

  A numeric vector of times at which to predict survival.

- type:

  Character string specifying the prediction output. Use `"event"` for
  the event survival matrix, `"censoring"` for the censoring survival
  matrix, or `"both"` for the full list of outputs.

- onlySL:

  Logical. If TRUE, only uses models with weights \> threshold.

- threshold:

  Numeric. The weight threshold for onlySL.

- ...:

  Additional ignored arguments.

## Value

If `type = "event"` or `type = "censoring"`, a numeric matrix with rows
corresponding to observations and columns corresponding to `new.times`.
If `type = "both"`, a list containing:

- `event.predict`: A numeric matrix of final event survival predictions.

- `event.library.predict`: A 3D numeric array of event learner
  predictions.

- `cens.predict`: A numeric matrix of final censoring survival
  predictions.

- `cens.library.predict`: A 3D numeric array of censoring learner
  predictions.

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:10, , drop = FALSE]
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

  pred_event <- predict(
    object = fit,
    newdata = newX,
    new.times = new.times,
    type = "event"
  )

  dim(pred_event)
}
#> [1] 10  6
```
