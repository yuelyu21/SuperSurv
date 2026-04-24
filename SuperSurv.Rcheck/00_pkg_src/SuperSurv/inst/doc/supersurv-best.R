## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
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

## ----fit-models, results='hide', message=FALSE, warning=FALSE-----------------
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

## ----inspect-weights----------------------------------------------------------
cat("\n--- ENSEMBLE WEIGHTS (selection = 'ensemble') ---\n")
print(round(event_weights(fit_ensemble), 4))

cat("\n--- BEST MODEL WEIGHTS (selection = 'best') ---\n")
print(round(event_weights(fit_best), 4))

## ----evaluate-comparison------------------------------------------------------
# Evaluate the Ensemble
cat("Performance of the ENSEMBLE Super Learner:\n")
eval_summary(fit_ensemble, newdata = X_te, time = test$duration, 
             event = test$event, eval_times = new.times)

cat("\nPerformance of the BEST MODEL Super Learner:\n")
eval_summary(fit_best, newdata = X_te, time = test$duration, 
             event = test$event, eval_times = new.times)

