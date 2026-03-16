# Keep All Variables Screener

Keep All Variables Screener

## Usage

``` r
screen.all(X, ...)
```

## Arguments

- X:

  Training covariate data.frame.

- ...:

  Additional ignored arguments.

## Value

A logical vector of the same length as the number of columns in `X`,
indicating which variables passed the screening algorithm (`TRUE` to
keep, `FALSE` to drop).

## Examples

``` r
data("metabric", package = "SuperSurv")
dat <- metabric[1:20, ]
x_cols <- grep("^x", names(dat))[1:5]
X <- dat[, x_cols, drop = FALSE]

screen.all(X)
#>   x0   x1   x2   x3   x4 
#> TRUE TRUE TRUE TRUE TRUE 
```
