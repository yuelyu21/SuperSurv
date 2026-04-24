# Access variables selected by SuperSurv screeners

Returns the variables retained by a screening step for one or more
fitted event or censoring learners.

## Usage

``` r
selected_variables(object, ...)

# S3 method for class 'SuperSurv'
selected_variables(object, type = c("event", "censoring"), learner = NULL, ...)
```

## Arguments

- object:

  A fitted object of class `"SuperSurv"`.

- ...:

  Additional arguments ignored.

- type:

  Character string specifying whether to inspect event or censoring
  learners.

- learner:

  Optional learner index or learner name. If omitted, selected variables
  are returned for every learner of the requested type.

## Value

A named list of character vectors, or a single character vector when
`learner` has length one.
