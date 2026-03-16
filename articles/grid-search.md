# 5. Advanced Hyperparameter Tuning & Grid Search

## Introduction

A single machine learning configuration is rarely optimal for every
dataset. A Random Forest with deep trees might overfit on small data,
while shallow trees might underfit on complex data.

In the Super Learner framework, we handle “tuning” by adding multiple
versions of the same algorithm to our library, each with different
hyperparameters. The meta-learner automatically evaluates them via
cross-validation and assigns the highest weight to the best
configuration.

`SuperSurv` provides the `create_grid` function to automate this process
entirely.

## 1. Prepare the Data

``` r
library(SuperSurv)
library(survival)

data("metabric", package = "SuperSurv")
set.seed(123)

train_idx <- sample(1:nrow(metabric), 0.7 * nrow(metabric))
train <- metabric[train_idx, ]
test  <- metabric[-train_idx, ]

X_tr <- train[, grep("^x", names(metabric))]
X_te <- test[, grep("^x", names(metabric))]

new.times <- seq(50, 200, by = 25) 
```

## 2. Defining the Search Space

Let’s assume we want to tune a Random Survival Forest (`surv.rfsrc`).
Instead of writing manual wrapper functions for every combination of
parameters, we define a list of the values we want to test.

``` r

# Define our hyperparameter search space
rf_tuning_params <- list(
  nodesize = c(3, 15, 30), # Deep vs. Shallow trees
  ntree = c(100, 500)      # Number of trees in the forest
)
```

## 3. Generating the Tuning Grid

We pass our base learner name and parameter list into
[`create_grid()`](https://yuelyu21.github.io/SuperSurv/reference/create_grid.md).
The function will automatically generate unique R functions for all 6
combinations in your global environment and return their names!

``` r
# Automatically generate 6 different Random Forest wrappers
rf_grid <- create_grid(
  base_learner = "surv.rfsrc", 
  grid_params = rf_tuning_params
)

# View the dynamically generated function names
print(rf_grid)
#> [1] "surv.rfsrc_nodesize3_ntree100"  "surv.rfsrc_nodesize15_ntree100"
#> [3] "surv.rfsrc_nodesize30_ntree100" "surv.rfsrc_nodesize3_ntree500" 
#> [5] "surv.rfsrc_nodesize15_ntree500" "surv.rfsrc_nodesize30_ntree500"
```

## 4. Fitting the Tuned Ensemble

You can now pass this dynamically generated character vector directly
into your `SuperSurv` library. We can combine our tuned ML grid with a
baseline Cox model to ensure we always have a linear fallback.

``` r
# Combine baseline with our 6 tuned forests
my_ml_library <- c("surv.coxph", rf_grid)

# Fit the ensemble (pseudo-code)
fit_tuned <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_ml_library,
  cens.library = c("surv.coxph"),
  nFolds = 3
)

# Inspect which hyperparameter combination won!
round(fit_tuned$event.coef, 3)
```

By leveraging
[`create_grid()`](https://yuelyu21.github.io/SuperSurv/reference/create_grid.md),
you can easily execute a massive hyperparameter search space across
XGBoost, SVMs, or Random Forests, mathematically guaranteeing that your
final ensemble minimizes the cross-validated risk without any manual
trial-and-error!
