## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----eval=FALSE---------------------------------------------------------------
# install.packages(c("future", "future.apply"))

## ----parallel-setup, message=FALSE, warning=FALSE, eval=FALSE-----------------
# library(SuperSurv)
# library(future)
# library(survival)
# 
# data("metabric", package = "SuperSurv")
# 
# # 1. Define the parallel plan
# # 'multisession' opens background R sessions.
# # We tell it to use 4 CPU cores (workers).
# plan(multisession, workers = 4)

## ----run-parallel, eval=FALSE-------------------------------------------------
# X <- metabric[, grep("^x", names(metabric))]
# new.times <- seq(50, 200, by = 25)
# 
# # 2. Run the model with parallel = TRUE
# fit_parallel <- SuperSurv(
#   time = metabric$duration,
#   event = metabric$event,
#   X = X,
#   newX = X,
#   new.times = new.times,
#   event.library = c("surv.coxph", "surv.weibull", "surv.rfsrc"),
#   cens.library = c("surv.coxph"),
#   parallel = TRUE,     # <--- The magic argument
#   nFolds = 5
# )

## ----close-parallel, eval=FALSE-----------------------------------------------
# # 3. Return to sequential processing
# plan(sequential)

