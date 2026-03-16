# Random Survival Forest Screening Algorithm

This screening algorithm uses the `randomForestSRC` package to select
covariates based on their Variable Importance (VIMP). It grows a fast
forest and retains features with a VIMP greater than zero.

## Usage

``` r
screen.rfsrc(
  time,
  event,
  X,
  obsWeights = NULL,
  minscreen = 2,
  ntree = 100,
  ...
)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- X:

  Training covariate data.frame or matrix.

- obsWeights:

  Numeric vector of observation weights.

- minscreen:

  Integer. Minimum number of covariates to return. Defaults to 2.

- ntree:

  Integer. Number of trees to grow. Defaults to 100 for fast screening.

- ...:

  Additional arguments passed to
  [`rfsrc`](https://www.randomforestsrc.org//reference/rfsrc.html).

## Value

A logical vector of the same length as the number of columns in `X`,
indicating which variables passed the screening algorithm (`TRUE` to
keep, `FALSE` to drop).

## Examples

``` r
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]

  screen.rfsrc(
    time = dat$duration,
    event = dat$event,
    X = X,
    minscreen = 2,
    ntree = 10
  )
}
#>    x0    x1    x2    x3    x4 
#>  TRUE  TRUE FALSE FALSE FALSE 
```
