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
set.seed(123)
train_idx <- sample(1:nrow(metabric), 0.7 * nrow(metabric))
train <- metabric[train_idx, ]
test  <- metabric[-train_idx, ]

X_tr <- train[, grep("^x", names(metabric))]
X_te <- test[, grep("^x", names(metabric))]

# Our max follow-up is well beyond 200, so this grid is safe.
new.times <- seq(50, 200, by = 25) 

## ----fit-models, results='hide', message=FALSE, warning=FALSE-----------------
my_library <- c("surv.coxph", "surv.weibull", "surv.rpart")

fit_supersurv <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_library,
  cens.library = my_library,
  control = list(saveFitLibrary = TRUE),
  verbose = FALSE,
  nFolds = 3
)

## ----evaluate-performance-----------------------------------------------------
# Evaluate performance directly using the fitted model and test data
performance_results <- eval_summary(
  object = fit_supersurv,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_times = new.times
)

## ----plot-benchmark, fig.align='center', fig.width=9, fig.height=3------------
plot_benchmark(
  object = fit_supersurv,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_times = new.times
)

## ----plot-calibration, fig.align='center'-------------------------------------
# Assess calibration specifically at Time = 150
plot_calibration(
  object = fit_supersurv,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_time = 150,
  bins = 5 # Group patients into 5 risk quintiles
)

