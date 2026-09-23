# Create a Tuning Grid of Survival Learners

Dynamically generates custom wrapper functions for a specified base
learner across a grid of hyperparameters.

## Usage

``` r
create_grid(base_learner, grid_params)
```

## Arguments

- base_learner:

  Character string of the base learner function name (e.g., "surv.gbm").

- grid_params:

  List of numeric/character vectors containing hyperparameter values.

## Value

A character vector of class `"SuperSurv_grid"` containing the newly
generated function names.

## Details

Generated functions are assigned to the calling environment. The base
learner is resolved there (with package functions as a fallback) when
the grid is created and retained by each generated function. Local base
learners and grids can therefore be used inside a function that calls
[`SuperSurv()`](https://yuelyu21.github.io/SuperSurv/reference/SuperSurv.md),
without modifying the global environment.
