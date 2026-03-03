# Calculate Restricted Mean Survival Time (RMST)

Calculate Restricted Mean Survival Time (RMST)

## Usage

``` r
get_rmst(surv_matrix, times, tau)
```

## Arguments

- surv_matrix:

  Matrix of survival probabilities (rows: patients, cols: time points)

- times:

  Vector of time points corresponding to the columns

- tau:

  The restriction time horizon

## Value

A vector of RMST values for each patient
