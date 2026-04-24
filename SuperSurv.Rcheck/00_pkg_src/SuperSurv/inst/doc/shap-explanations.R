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

my_library <- c("surv.coxph", "surv.weibull", "surv.rfsrc")

fit_sl <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_library,
  cens.library = c("surv.coxph"),
  verbose = FALSE,
  selection = "ensemble",
  nFolds = 3
)

## ----calc-shap, message=FALSE, warning=FALSE----------------------------------
# Explain the first 50 test patients using 100 training patients as the background
X_explain_subset <- X_te[1:50, ]
X_background_subset <- X_tr[1:100, ]

# Calculate weighted Ensemble SHAP values using Kernel SHAP
shap_vals <- explain_kernel(
  model = fit_sl, 
  X_explain = X_explain_subset, 
  X_background = X_background_subset, 
  nsim = 20
)

## ----plot-global, fig.width=6, fig.height=4-----------------------------------
plot_global_importance(shap_vals, top_n = 5)

## ----plot-beeswarm, fig.width=7, fig.height=5, warning=FALSE------------------
plot_beeswarm(shap_vals, data = X_explain_subset, top_n = 5)

## ----plot-waterfall, fig.width=6, fig.height=4--------------------------------
# Explain Patient #1 from our test subset
plot_patient_waterfall(shap_vals, patient_index = 1, top_n = 5)

## ----plot-dependence, fig.width=6, fig.height=4, warning=FALSE----------------
plot_dependence(shap_vals, data = X_explain_subset, feature_name = "x0")

## ----plot-heatmap, fig.width=7, fig.height=5----------------------------------
# Plot the survival trajectories for the first 50 test patients
plot_survival_heatmap(fit_sl, newdata = X_explain_subset, times = new.times)

## ----survex-bridge, message=FALSE, warning=FALSE------------------------------
library(survex)

# 1. Create the true survival object for the explanation subset
y_explain <- survival::Surv(test$duration[1:50], test$event[1:50])

# 2. Build the survex explainer using our custom function
surv_explainer <- explain_survex(
  model = fit_sl, 
  data = X_explain_subset, 
  y = y_explain, 
  times = new.times
)

## ----survex-importance, fig.width=7, fig.height=5, message=FALSE, warning=FALSE----
# Calculate time-dependent model parts (permutation feature importance)
time_importance <- model_parts(surv_explainer)

# Plot the dynamic importance over time
plot(time_importance)

## ----survex-profile, fig.width=7, fig.height=5, message=FALSE, warning=FALSE----
# Calculate the partial dependence profile for feature 'x0'
pdp_time <- model_profile(surv_explainer, variables = "x0")

plot(pdp_time)

## ----survex-survshap, fig.width=7, fig.height=5, message=FALSE, warning=FALSE----
# Explain Patient #1 over time
patient_1_data <- X_explain_subset[1, , drop = FALSE]

# Calculate SurvSHAP(t)
survshap_t <- predict_parts(surv_explainer, new_observation = patient_1_data, type = "survshap")

plot(survshap_t)

## ----survex-performance, fig.width=7, fig.height=5, message=FALSE, warning=FALSE----
# Calculate time-dependent performance metrics via survex
survex_perf <- model_performance(surv_explainer)

# Plot the Brier score and AUC curves
plot(survex_perf)

