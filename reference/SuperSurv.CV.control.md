# Control parameters for Cross-Validation in SuperSurv

Control parameters for Cross-Validation in SuperSurv

## Usage

``` r
SuperSurv.CV.control(
  V = 10L,
  stratifyCV = TRUE,
  shuffle = TRUE,
  validRows = NULL
)
```

## Arguments

- V:

  Number of folds. Default is 10.

- stratifyCV:

  Logical. If TRUE, ensures event rates are balanced across folds.

- shuffle:

  Logical. If TRUE, shuffles rows before splitting.

- validRows:

  Optional custom list of indices for folds.

## Value

A list of CV parameters.
