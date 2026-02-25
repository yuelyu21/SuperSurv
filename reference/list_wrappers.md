# List Available Wrappers and Screeners in SuperSurv

This function prints all built-in prediction algorithms and feature
screening algorithms available in the `SuperSurv` package.

## Usage

``` r
list_wrappers(what = "both")
```

## Arguments

- what:

  Character string. If `"both"` (default), lists both prediction and
  screening functions. If `"surv"`, lists only prediction models. If
  `"screen"`, lists only screening algorithms. Otherwise, lists all
  exports.

## Value

An invisible character vector containing the requested function names.
