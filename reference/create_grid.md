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

A character vector of the newly generated function names.
