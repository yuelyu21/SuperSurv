# Access SuperSurv ensemble weights

Extracts the fitted event or censoring ensemble weights from a
`SuperSurv` object.

## Usage

``` r
event_weights(object, ...)

# S3 method for class 'SuperSurv'
event_weights(object, ...)

censor_weights(object, ...)

# S3 method for class 'SuperSurv'
censor_weights(object, ...)
```

## Arguments

- object:

  A fitted object of class `"SuperSurv"`.

- ...:

  Additional arguments ignored.

## Value

A named numeric vector of ensemble weights.

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

  event_weights(fit)
  censor_weights(fit)
}
#> surv.coxph_screen.all 
#>                     1 
```
