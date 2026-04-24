# Extending SuperSurv

This vignette describes the practical contract used by
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md)
to call survival learners and screening methods. The goal is to make
small extensions possible without changing the ensemble algorithm.

## How Learners Are Found

Learners are supplied to
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md)
by name:

``` r
event.library <- c("surv.coxph", "my.surv.km")
```

During cross-validation and final fitting,
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md)
looks up each name as an R function and calls it directly. A custom
learner can therefore live in a user script, package, or analysis file,
as long as it is available in the R session before calling
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md).

## Learner Wrapper Contract

A survival learner wrapper should have this shape:

``` r
my.surv.learner <- function(time, event, X, newdata, new.times,
                            obsWeights, id, ...) {
  # fit model on X, time, event
  # predict survival probabilities for newdata at new.times
  list(pred = pred_matrix, fit = fitted_object)
}
```

The required arguments are:

| Argument     | Meaning                                                       |
|--------------|---------------------------------------------------------------|
| `time`       | numeric observed follow-up time                               |
| `event`      | event indicator, with `1` for the event and `0` for censoring |
| `X`          | training covariates as a data frame                           |
| `newdata`    | covariates for prediction as a data frame                     |
| `new.times`  | time grid where survival probabilities are needed             |
| `obsWeights` | observation weights, possibly `NULL`                          |
| `id`         | optional cluster or subject identifier, possibly `NULL`       |
| `...`        | learner-specific tuning parameters                            |

The returned `pred` must be a numeric matrix with:

- one row per row of `newdata`,
- one column per value of `new.times`,
- entries on the survival-probability scale, usually between `0` and
  `1`.

If possible, each row should be non-increasing across time. Wrappers
commonly clamp small numerical drift back into `[0, 1]`.

The returned `fit` can be any object needed for future prediction. When
`control = list(saveFitLibrary = TRUE)`,
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md)
stores these fitted objects and
[`predict.SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/predict.SuperSurv.md)
later calls [`predict()`](https://rdrr.io/r/stats/predict.html) on them.

## Prediction Method Contract

If a custom learner should support prediction on new data after fitting,
store a classed object in `fit` and define a matching S3 prediction
method:

``` r
predict.my.surv.learner <- function(object, newdata, new.times, ...) {
  pred_matrix
}
```

This method should return the same type of matrix as the wrapper’s
`pred` component: `nrow(newdata)` rows by `length(new.times)` columns.

## Direct Survival Curves vs Risk Scores

Some survival learners directly return survival curves. In that case,
the wrapper usually interpolates the native prediction grid to
`new.times`.

Other learners return a risk score or linear predictor. Those wrappers
need one extra calibration step:

1.  fit the risk-score model,
2.  estimate a baseline cumulative hazard on the training data,
3.  predict risk scores for `newdata`,
4.  convert scores to survival curves with
    `S(t | x) = exp(-exp(eta(x)) H0(t))`.

Several built-in Cox-style wrappers follow this pattern. The internal
helper `safe_breslow_step()` is used by some package wrappers, but it is
not exported; custom wrappers should either use a standard package
prediction method that already returns survival curves or implement
their own baseline-hazard calibration carefully.

## Minimal Custom Learner

The following learner ignores covariates and fits a weighted
Kaplan-Meier curve. It is intentionally simple, but it satisfies the
full
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md)
learner contract.

``` r
my.surv.km <- function(time, event, X, newdata, new.times,
                       obsWeights = NULL, id = NULL, ...) {
  if (is.null(obsWeights)) {
    obsWeights <- rep(1, length(time))
  }

  fit_km <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    weights = obsWeights
  )

  step_surv <- stats::stepfun(fit_km$time, c(1, fit_km$surv), right = FALSE)
  surv_probs <- step_surv(new.times)

  pred <- matrix(
    surv_probs,
    nrow = nrow(newdata),
    ncol = length(new.times),
    byrow = TRUE
  )

  fit <- list(object = fit_km)
  class(fit) <- "my.surv.km"

  list(pred = pred, fit = fit)
}

predict.my.surv.km <- function(object, newdata, new.times, ...) {
  fit_km <- object$object
  step_surv <- stats::stepfun(fit_km$time, c(1, fit_km$surv), right = FALSE)
  surv_probs <- step_surv(new.times)

  matrix(
    surv_probs,
    nrow = nrow(newdata),
    ncol = length(new.times),
    byrow = TRUE
  )
}
```

You can test the wrapper outside the ensemble before adding it to a
library:

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
times <- seq(50, 150, by = 50)

wrapper_out <- my.surv.km(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X[1:5, , drop = FALSE],
  new.times = times
)

dim(wrapper_out[["pred"]])
#> [1] 5 3
dim(predict(wrapper_out[["fit"]], newdata = X[1:5, , drop = FALSE], new.times = times))
#> [1] 5 3
```

The custom learner can be mixed with built-in learners:

``` r
fit <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_train,
  newdata = X_test,
  new.times = seq(50, 300, by = 50),
  event.library = c("surv.coxph", "my.surv.km"),
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE)
)

event_weights(fit)
```

## Screening Method Contract

Screening methods are also supplied by name. A screener should accept
the same basic training inputs and return a logical vector with one
value per column of `X`:

``` r
my.screen.first3 <- function(time, event, X, obsWeights = NULL, id = NULL, ...) {
  keep <- rep(FALSE, ncol(X))
  keep[seq_len(min(3, ncol(X)))] <- TRUE
  names(keep) <- names(X)
  keep
}
```

Use a list entry to pair a learner with one or more screeners:

``` r
fit <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_train,
  newdata = X_test,
  new.times = seq(50, 300, by = 50),
  event.library = list(c("surv.coxph", "screen.all", "my.screen.first3")),
  cens.library = c("surv.coxph")
)

selected_variables(fit, learner = 2)
```

## Development Checklist

Before using a custom wrapper in a large ensemble, check:

- `pred` is a numeric matrix with dimensions `nrow(newdata)` by
  `length(new.times)`,
- predictions are on the survival-probability scale,
- `predict(wrapper_out[["fit"]], newdata, new.times)` returns the same
  shape when saved fits are needed,
- the wrapper tolerates `obsWeights = NULL` and `id = NULL`,
- all required packages are checked with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html),
- any screener returns a logical vector aligned with the columns of `X`.

These rules are deliberately small. They allow new learners to be added
without changing the SuperSurv ensemble machinery.
