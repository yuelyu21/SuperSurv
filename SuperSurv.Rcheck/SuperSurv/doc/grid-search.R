## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
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

## ----message=FALSE, warning=FALSE---------------------------------------------

# Define our hyperparameter search space
rf_tuning_params <- list(
  nodesize = c(3, 15, 30), # Deep vs. Shallow trees
  ntree = c(100, 500)      # Number of trees in the forest
)

## ----create-grid--------------------------------------------------------------
# Automatically generate 6 different Random Forest wrappers
rf_grid <- create_grid(
  base_learner = "surv.rfsrc", 
  grid_params = rf_tuning_params
)

# View the dynamically generated function names
print(rf_grid)

## ----fit-tuned-model, eval=FALSE----------------------------------------------
# # Combine baseline with our 6 tuned forests
# my_ml_library <- c("surv.coxph", rf_grid)
# 
# # Fit the ensemble (pseudo-code)
# fit_tuned <- SuperSurv(
#   time = train$duration,
#   event = train$event,
#   X = X_tr,
#   newdata = X_te,
#   new.times = new.times,
#   event.library = my_ml_library,
#   cens.library = c("surv.coxph"),
#   nFolds = 3
# )
# 
# # Inspect which hyperparameter combination won!
# round(event_weights(fit_tuned), 3)

