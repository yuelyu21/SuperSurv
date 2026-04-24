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

train_idx <- sample(1:nrow(metabric), 0.7 * nrow(metabric))
train <- metabric[train_idx, ]
test  <- metabric[-train_idx, ]

X_tr <- train[, grep("^x", names(metabric))]
X_te <- test[, grep("^x", names(metabric))]
new.times <- seq(50, 200, by = 25)

# Define a library covering different parametric assumptions
parametric_library <- c("surv.coxph",       # Semi-parametric baseline
                        "surv.weibull",     # Assumes hazard increases/decreases monotonically
                        "surv.exponential", # Assumes constant hazard over time
                        "surv.lognormal")   # Assumes hazard rises then falls

## ----fit-parametric, results='hide', message=FALSE, warning=FALSE-------------
fit_parametric <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = parametric_library,
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE),
  verbose = FALSE,
  selection = "ensemble",
  nFolds = 3
)

## ----evaluate-parametric------------------------------------------------------
summary(fit_parametric)

