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

new.times <- seq(50, 200, by = 25) 

## ----standalone-rf------------------------------------------------------------
# 1. Fit the standalone wrapper
rf_standalone <- surv.rfsrc(
  time = train$duration,
  event = train$event,
  X = X_tr,
  new.times = new.times 
)

# 2. Extract the fitted model object and prediction matrix
rf_fit <- rf_standalone[["fit"]]
rf_pred_matrix <- rf_standalone[["pred"]]

## ----plot-standalone, fig.align='center'--------------------------------------
# Plot the first 3 patients in our training set
plot_predict(preds = rf_pred_matrix, eval_times = new.times, patient_idx = 1:3)

## ----eval-standalone, fig.align='center', fig.height=4, fig.width= 9----------
# The function automatically detects this is a single model and plots it!
plot_benchmark(
  object = rf_fit,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_times = new.times
)

## -----------------------------------------------------------------------------

# plot_calibration(
#   object   = rf_fit,
#   newdata  = X_te,
#   time     = test$duration,
#   event    = test$event,
#   eval_time = 150,
#   bins     = 2
# )


## ----fit-models, results='hide', message=FALSE, warning=FALSE-----------------
my_library <- c("surv.coxph", "surv.weibull", "surv.rfsrc")

fit_supersurv <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_library,
  cens.library = c("surv.coxph"), 
  control = list(saveFitLibrary = TRUE),
  verbose = FALSE,
  nFolds = 3
)

## ----plot-ensemble-benchmark, fig.align='center', fig.height=9----------------
plot_benchmark(
  object = fit_supersurv,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_times = new.times
)

