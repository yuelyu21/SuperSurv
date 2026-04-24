# Summarize a SuperSurv fit

Summarizes the fitted event and censoring ensembles, including learner
names, ensemble weights, cross-validated risks, error flags, prediction
dimensions, and recorded timing information.

## Usage

``` r
# S3 method for class 'SuperSurv'
summary(object, ...)

# S3 method for class 'summary.SuperSurv'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  A fitted object of class `"SuperSurv"`.

- ...:

  Additional arguments ignored.

- x:

  A summary object produced by `summary.SuperSurv()`.

- digits:

  Number of significant digits to use for displayed weights and risks.

## Value

An object of class `"summary.SuperSurv"`, a list containing the matched
call, selection mode, event and censoring learner summaries, prediction
dimensions, and timing information.
