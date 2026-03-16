# Marginal Cox Regression Screening

Marginal Cox Regression Screening

## Usage

``` r
screen.marg(time, event, X, obsWeights = NULL, minscreen = 2, min.p = 0.1, ...)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data.frame.

- obsWeights:

  Observation weights.

- minscreen:

  Minimum number of covariates to return. Defaults to 2.

- min.p:

  Threshold p-value. Defaults to 0.1.

- ...:

  Additional ignored arguments.

## Value

A logical vector of the same length as the number of columns in `X`,
indicating which variables passed the screening algorithm (`TRUE` to
keep, `FALSE` to drop).

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:5]
X <- dat[, x_cols, drop = FALSE]

screen.marg(
  time = dat$duration,
  event = dat$event,
  X = X,
  minscreen = 2,
  min.p = 0.2
)
#>    x0    x1    x2    x3    x4 
#>  TRUE  TRUE FALSE FALSE FALSE 
```
