# Kaplan-Meier Prediction Algorithm

This prediction algorithm ignores all covariates and computes the
marginal Kaplan-Meier survival estimator using the
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) function.

## Usage

``` r
surv.km(time, event, X, newdata, new.times, obsWeights, id, ...)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- X:

  Training covariate data.frame (Ignored by KM).

- newdata:

  Test covariate data.frame to use for prediction.

- new.times:

  Numeric vector of times at which to predict survival.

- obsWeights:

  Numeric vector of observation weights.

- id:

  Optional vector indicating subject/cluster identities.

- ...:

  Additional ignored arguments.

## Value

A list containing:

- `fit`: A list containing the fitted
  [`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) object.

- `pred`: A numeric matrix of cross-validated survival predictions
  evaluated at `new.times`.

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.km(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit$pred)
#> [1] 5 3
```
