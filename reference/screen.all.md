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
