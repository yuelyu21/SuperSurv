# 3. Ensemble vs. Best Model Selection

## Introduction

In the theoretical framework of Super Learning, there are two distinct
ways to utilize the cross-validated risks of your base algorithms:

1.  **The Ensemble Super Learner:** Calculates a weighted average
    (convex combination) of the base learners. This “soft selection”
    smooths out variance, integrates different feature spaces, and
    generally yields the lowest finite-sample prediction error.
2.  **The “Best” Model Selector:** Identifies the single algorithm with
    the lowest cross-validated risk and assigns it a weight of `1.0`
    (and all others `0.0`). This is a “hard selection” or
    “winner-take-all” approach.

While the Ensemble is asymptotically optimal, selecting the single best
model is incredibly useful when interpretability is strictly tied to one
specific algorithm family. Instead of manually cherry-picking the best
model—which introduces researcher bias and invalidates post-selection
inference—`SuperSurv` automates the selection using rigorous, internal
cross-validation.

You can easily toggle between these two paradigms using the `selection`
argument.

## 1. Prepare the Data

We load the `metabric` dataset and define our evaluation time grid
exactly as we did in the previous tutorials.

``` r
library(SuperSurv)
library(survival)

data("metabric", package = "SuperSurv")
set.seed(42)

# Standard Train/Test Split
train_idx <- sample(1:nrow(metabric), 0.7 * nrow(metabric))
train <- metabric[train_idx, ]
test  <- metabric[-train_idx, ]

X_tr <- train[, grep("^x", names(metabric))]
X_te <- test[, grep("^x", names(metabric))]
new.times <- seq(50, 200, by = 25)

# Define a diverse library of base learners
my_library <- c("surv.coxph", "surv.weibull", "surv.rpart")
```

## 2. Fit Both Super Learners

We will fit two separate `SuperSurv` models. The first will use the
default Ensemble approach (`selection = "ensemble"`), and the second
will use the winner-take-all approach (`selection = "best"`).

``` r
# 1. The Ensemble Super Learner (Weighted Average)
fit_ensemble <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_library,
  cens.library = my_library,
  control = list(saveFitLibrary = TRUE),
  verbose = FALSE,
  selection = "ensemble",     # <-- Calculates fractional weights
  nFolds = 3
)

# 2. The 'Best' Super Learner (Winner-Take-All)
fit_best <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_library,
  cens.library = my_library,
  control = list(saveFitLibrary = TRUE),
  verbose = FALSE,
  selection = "best",         # <-- Selects the single best model
  nFolds = 3
)
```

## 3. Inspecting the Weights

The difference between the two methodologies is immediately obvious when
we inspect the meta-learner coefficients (`event.coef`).

``` r
cat("\n--- ENSEMBLE WEIGHTS (selection = 'ensemble') ---\n")
#> 
#> --- ENSEMBLE WEIGHTS (selection = 'ensemble') ---
print(round(fit_ensemble$event.coef, 4))
#>   surv.coxph_screen.all surv.weibull_screen.all   surv.rpart_screen.all 
#>                  0.5152                  0.1680                  0.3167

cat("\n--- BEST MODEL WEIGHTS (selection = 'best') ---\n")
#> 
#> --- BEST MODEL WEIGHTS (selection = 'best') ---
print(round(fit_best$event.coef, 4))
#> [1] 1 0 0
```

**Interpretation:** Notice that the Ensemble Super Learner distributes
the weight across multiple models to minimize the overall loss function.
The “Best” Super Learner simply looks at the cross-validated risk,
identifies the champion, and gives it 100% of the weight.

## 4. Evaluate and Compare Performance

A common question in clinical research is: *“Does the complexity of the
Ensemble actually perform better than simply picking the best single
model?”*

We can answer this instantly using
[`eval_summary()`](https://yuelyu21.github.io/SuperSurv/reference/eval_summary.md).

``` r
# Evaluate the Ensemble
cat("Performance of the ENSEMBLE Super Learner:\n")
#> Performance of the ENSEMBLE Super Learner:
eval_summary(fit_ensemble, newdata = X_te, time = test$duration, 
             event = test$event, eval_times = new.times)
#>                     Model    IBS  Uno_C   iAUC
#> 1      SuperSurv_Ensemble 0.1978 0.6569 0.6788
#> 2   surv.coxph_screen.all 0.1989 0.6525 0.6744
#> 3 surv.weibull_screen.all 0.1993 0.6532 0.6755
#> 4   surv.rpart_screen.all 0.2068 0.6083 0.6541

cat("\nPerformance of the BEST MODEL Super Learner:\n")
#> 
#> Performance of the BEST MODEL Super Learner:
eval_summary(fit_best, newdata = X_te, time = test$duration, 
             event = test$event, eval_times = new.times)
#>                     Model    IBS  Uno_C   iAUC
#> 1      SuperSurv_Ensemble 0.1989 0.6525 0.6744
#> 2   surv.coxph_screen.all 0.1989 0.6525 0.6744
#> 3 surv.weibull_screen.all 0.1993 0.6532 0.6755
#> 4   surv.rpart_screen.all 0.2068 0.6083 0.6541
```

By comparing the resulting Brier scores and C-indices, you can
empirically justify whether the “soft selection” of the ensemble is
mathematically necessary for your specific dataset, or if the “hard
selection” of a single model is sufficient.
