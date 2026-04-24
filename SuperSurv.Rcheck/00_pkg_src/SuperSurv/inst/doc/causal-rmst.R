## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(SuperSurv)
set.seed(123)

# Load built-in data
data("metabric", package = "SuperSurv")

# Define predictors and time grid
X <- metabric[, grep("^x", names(metabric))]
new.times <- seq(10, 150, by = 10)

## ----train-model--------------------------------------------------------------
fit <- SuperSurv(
  time = metabric$duration,
  event = metabric$event,
  X = X,
  newdata = X,
  new.times = new.times,
  event.library = c("surv.coxph", "surv.rfsrc"),
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE) 
)

## ----causal-effect------------------------------------------------------------
# Estimate the adjusted difference up to tau = 100 months
results <- estimate_marginal_rmst(
  fit = fit, 
  data = metabric, 
  trt_col = "x4", 
  times = new.times, 
  tau = 100
)



print(results$ATE_RMST)

## -----------------------------------------------------------------------------
rmst_results_inf <- estimate_marginal_rmst(
  fit = fit,
  data = metabric,
  trt_col = "x4",
  times = new.times,
  tau = 100,
  inference = TRUE,
  B = 100,
  seed = 123
)

rmst_results_inf$ATE_RMST
rmst_results_inf$SE_RMST
rmst_results_inf$CI_RMST
format.pval(rmst_results_inf$p_value, digits = 3, eps = 1e-16)

## ----plot-curve---------------------------------------------------------------
# Plot the Delta RMST across a sequence of tau values
tau_grid <- seq(20, 140, by = 30)
plot_marginal_rmst_curve(
  fit = fit, 
  data = metabric, 
  trt_col = "x4", 
  times = new.times, 
  tau_seq = tau_grid,
  inference = TRUE, 
  B = 100, 
  seed = 123, 
  ci_level = 0.95
)

## ----plot-obs-----------------------------------------------------------------
plot_rmst_vs_obs(
  fit = fit, 
  data = metabric, 
  time_col = "duration", 
  event_col = "event", 
  times = new.times, 
  tau = 350
)

