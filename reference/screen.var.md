# High Variance Screening Algorithm (Unsupervised)

An unsupervised screening algorithm that filters out low-variance
features. This is particularly useful for high-dimensional genomic or
transcriptomic data where many features remain relatively constant
across all observations.

## Usage

``` r
screen.var(
  time,
  event,
  X,
  obsWeights = NULL,
  keep_fraction = 0.5,
  minscreen = 2,
  ...
)
```

## Arguments

- time:

  Numeric vector of observed follow-up times (Ignored internally).

- event:

  Numeric vector of event indicators (Ignored internally).

- X:

  Training covariate data.frame or matrix.

- obsWeights:

  Numeric vector of observation weights (Ignored internally).

- keep_fraction:

  Numeric value between 0 and 1. The fraction of highest-variance
  features to retain. Defaults to 0.5 (keeps the top 50%).

- minscreen:

  Integer. Minimum number of covariates to return. Defaults to 2.

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
x_cols <- grep("^x", names(dat))[1:6]
X <- dat[, x_cols, drop = FALSE]

screen.var(
  time = dat$duration,
  event = dat$event,
  X = X,
  keep_fraction = 0.5,
  minscreen = 2
)
#>    x0    x1    x2    x3    x4    x5 
#>  TRUE  TRUE  TRUE FALSE FALSE FALSE 
```
