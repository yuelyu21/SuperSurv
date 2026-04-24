# Extract SuperSurv ensemble coefficients

Returns the optimized convex-combination weights from a fitted
`SuperSurv` object.

## Usage

``` r
# S3 method for class 'SuperSurv'
coef(object, type = c("event", "censoring", "both"), ...)
```

## Arguments

- object:

  A fitted object of class `"SuperSurv"`.

- type:

  Character string specifying which weights to return. Use `"event"` for
  the event-time ensemble, `"censoring"` for the censoring ensemble, or
  `"both"` for both sets of weights.

- ...:

  Additional arguments ignored.

## Value

For `type = "event"` or `type = "censoring"`, a named numeric vector of
ensemble weights. For `type = "both"`, a list with elements `event` and
`censoring`.

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
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

  coef(fit)
  coef(fit, type = "both")
}
#> $event
#> surv.coxph_screen.all surv.ridge_screen.all 
#>             0.1686131             0.8313869 
#> 
#> $censoring
#> surv.coxph_screen.all 
#>                     1 
#> 
```
