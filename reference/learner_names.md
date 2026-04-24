# Access SuperSurv learner names

Returns the fitted learner names from a `SuperSurv` object.

## Usage

``` r
learner_names(object, ...)

# S3 method for class 'SuperSurv'
learner_names(object, type = c("event", "censoring", "both"), ...)
```

## Arguments

- object:

  A fitted object of class `"SuperSurv"`.

- ...:

  Additional arguments ignored.

- type:

  Character string specifying whether to return event learner names,
  censoring learner names, or both.

## Value

For `type = "event"` or `type = "censoring"`, a character vector of
learner names. For `type = "both"`, a list with elements `event` and
`censoring`.
