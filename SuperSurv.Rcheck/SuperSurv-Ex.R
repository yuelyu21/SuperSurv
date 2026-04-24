pkgname <- "SuperSurv"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SuperSurv')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SuperSurv")
### * SuperSurv

flush(stderr()); flush(stdout())

### Name: SuperSurv
### Title: Super Learner for conditional survival functions
### Aliases: SuperSurv

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  fit
  event_weights(fit)
}



cleanEx()
nameEx("coef.SuperSurv")
### * coef.SuperSurv

flush(stderr()); flush(stdout())

### Name: coef.SuperSurv
### Title: Extract SuperSurv ensemble coefficients
### Aliases: coef.SuperSurv

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  coef(fit)
  coef(fit, type = "both")
}



cleanEx()
nameEx("estimate_marginal_rmst")
### * estimate_marginal_rmst

flush(stderr()); flush(stdout())

### Name: estimate_marginal_rmst
### Title: Estimate an Adjusted Marginal RMST Contrast
### Aliases: estimate_marginal_rmst

### ** Examples

## Not run: 
##D data("metabric", package = "SuperSurv")
##D x_cols <- grep("^x", names(metabric), value = TRUE)
##D X <- metabric[, x_cols]
##D new.times <- seq(10, 150, by = 10)
##D 
##D fit <- SuperSurv(
##D   time = metabric$duration,
##D   event = metabric$event,
##D   X = X,
##D   newdata = X,
##D   new.times = new.times,
##D   event.library = c("surv.coxph", "surv.rfsrc"),
##D   cens.library = c("surv.coxph"),
##D   control = list(saveFitLibrary = TRUE),
##D   nFolds = 3
##D )
##D 
##D rmst_res <- estimate_marginal_rmst(
##D   fit = fit,
##D   data = metabric,
##D   trt_col = "x4",
##D   times = new.times,
##D   tau = 100,
##D   inference = TRUE,
##D   B = 200,
##D   seed = 123
##D )
##D 
##D rmst_res$ATE_RMST
##D rmst_res$SE_RMST
##D rmst_res$CI_RMST
##D format.pval(rmst_res$p_value, digits = 3, eps = 1e-16)
## End(Not run)




cleanEx()
nameEx("eval_brier")
### * eval_brier

flush(stderr()); flush(stdout())

### Name: eval_brier
### Title: IPCW Brier Score and Integrated Brier Score (IBS)
### Aliases: eval_brier

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:10, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.coxph(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

eval_brier(
  time = dat$duration[1:10],
  event = dat$event[1:10],
  S_mat = fit[["pred"]],
  times = times
)



cleanEx()
nameEx("eval_cindex")
### * eval_cindex

flush(stderr()); flush(stdout())

### Name: eval_cindex
### Title: Calculate Concordance Index (Harrell's or Uno's)
### Aliases: eval_cindex

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:10, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.coxph(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

eval_cindex(
  time = dat$duration[1:10],
  event = dat$event[1:10],
  S_mat = fit[["pred"]],
  times = times,
  eval_time = 100,
  method = "uno"
)



cleanEx()
nameEx("eval_summary")
### * eval_summary

flush(stderr()); flush(stdout())

### Name: eval_summary
### Title: Evaluate SuperSurv predictions on test data
### Aliases: eval_summary

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:80, ]
x_cols <- grep("^x", names(dat))[1:3]

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = dat[, x_cols, drop = FALSE],
  new.times = seq(50, 200, by = 50),
  event.library = c("surv.coxph", "surv.km"),
  cens.library = c("surv.coxph", "surv.km")
)

res <- eval_summary(
  object = fit,
  newdata = dat[, x_cols, drop = FALSE],
  time = dat$duration,
  event = dat$event,
  eval_times = seq(50, 200, by = 50)
)

res




cleanEx()
nameEx("eval_timeROC")
### * eval_timeROC

flush(stderr()); flush(stdout())

### Name: eval_timeROC
### Title: Time-Dependent AUC and Integrated AUC
### Aliases: eval_timeROC

### ** Examples

 data("metabric", package = "SuperSurv")
 dat <- metabric[1:40, ]
 x_cols <- grep("^x", names(dat))[1:3]
 X <- dat[, x_cols, drop = FALSE]
 newX <- X[1:10, , drop = FALSE]
 times <- seq(50, 150, by = 50)

 fit <- surv.coxph(
   time = dat$duration,
   event = dat$event,
   X = X,
   newdata = newX,
   new.times = times,
   obsWeights = rep(1, nrow(dat)),
   id = NULL
 )

 eval_timeROC(
   time = dat$duration[1:10],
   event = dat$event[1:10],
   S_mat = fit[["pred"]],
   times = times
 )



cleanEx()
nameEx("event_weights")
### * event_weights

flush(stderr()); flush(stdout())

### Name: event_weights
### Title: Access SuperSurv ensemble weights
### Aliases: event_weights event_weights.SuperSurv censor_weights
###   censor_weights.SuperSurv

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  event_weights(fit)
  censor_weights(fit)
}



cleanEx()
nameEx("explain_kernel")
### * explain_kernel

flush(stderr()); flush(stdout())

### Name: explain_kernel
### Title: Explain Predictions with Global SHAP (Kernel SHAP)
### Aliases: explain_kernel

### ** Examples

if (requireNamespace("fastshap", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  shap_values <- explain_kernel(
    model = fit,
    X_explain = X[1:10, , drop = FALSE],
    X_background = X[11:40, , drop = FALSE],
    nsim = 5
  )

  dim(shap_values)
}



cleanEx()
nameEx("explain_survex")
### * explain_survex

flush(stderr()); flush(stdout())

### Name: explain_survex
### Title: Create a Time-Dependent Survex Explainer
### Aliases: explain_survex

### ** Examples

if (requireNamespace("survex", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  times <- seq(20, 120, by = 20)
  y <- survival::Surv(dat$duration, dat$event)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  explainer <- explain_survex(
    model = fit,
    data = X,
    y = y,
    times = times,
    label = "SuperSurv_demo"
  )

  class(explainer)
}



cleanEx()
nameEx("list_wrappers")
### * list_wrappers

flush(stderr()); flush(stdout())

### Name: list_wrappers
### Title: List Available Wrappers and Screeners in SuperSurv
### Aliases: list_wrappers

### ** Examples

list_wrappers()



cleanEx()
nameEx("plot_beeswarm")
### * plot_beeswarm

flush(stderr()); flush(stdout())

### Name: plot_beeswarm
### Title: Beeswarm Summary Plot for SuperSurv SHAP
### Aliases: plot_beeswarm

### ** Examples

if (requireNamespace("fastshap", quietly = TRUE) &&
    requireNamespace("ggforce", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  shap_values <- explain_kernel(
    model = fit,
    X_explain = X[1:20, , drop = FALSE],
    X_background = X[21:50, , drop = FALSE],
    nsim = 5
  )

  plot_beeswarm(
    shap_values = shap_values,
    data = X[1:20, , drop = FALSE],
    top_n = 5
  )
}



cleanEx()
nameEx("plot_benchmark")
### * plot_benchmark

flush(stderr()); flush(stdout())

### Name: plot_benchmark
### Title: Plot Longitudinal Benchmark Metrics
### Aliases: plot_benchmark

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:120, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  eval_times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = eval_times,
    event.library = c("surv.coxph", "surv.ranger"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  plot_benchmark(
    object = fit,
    newdata = X,
    time = dat$duration,
    event = dat$event,
    eval_times = eval_times,
    metrics = c("brier")
  )
}



cleanEx()
nameEx("plot_calibration")
### * plot_calibration

flush(stderr()); flush(stdout())

### Name: plot_calibration
### Title: Plot Survival Calibration Curve
### Aliases: plot_calibration

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:120, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  eval_times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = eval_times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  plot_calibration(
    object = fit,
    newdata = X,
    time = dat$duration,
    event = dat$event,
    eval_time = 100,
    bins = 4
  )
}



cleanEx()
nameEx("plot_dependence")
### * plot_dependence

flush(stderr()); flush(stdout())

### Name: plot_dependence
### Title: Plot SHAP Dependence for SuperSurv
### Aliases: plot_dependence

### ** Examples

if (requireNamespace("fastshap", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  shap_values <- explain_kernel(
    model = fit,
    X_explain = X[1:20, , drop = FALSE],
    X_background = X[21:50, , drop = FALSE],
    nsim = 5
  )

  plot_dependence(
    shap_values = shap_values,
    data = X[1:20, , drop = FALSE],
    feature_name = colnames(X)[1]
  )
}



cleanEx()
nameEx("plot_global_importance")
### * plot_global_importance

flush(stderr()); flush(stdout())

### Name: plot_global_importance
### Title: Plot Global Feature Importance for SuperSurv
### Aliases: plot_global_importance

### ** Examples

if (requireNamespace("fastshap", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  shap_values <- explain_kernel(
    model = fit,
    X_explain = X[1:10, , drop = FALSE],
    X_background = X[11:40, , drop = FALSE],
    nsim = 5
  )

  plot_global_importance(shap_values, top_n = 5)
}



cleanEx()
nameEx("plot_marginal_rmst_curve")
### * plot_marginal_rmst_curve

flush(stderr()); flush(stdout())

### Name: plot_marginal_rmst_curve
### Title: Plot Adjusted Marginal RMST Contrast Over Time
### Aliases: plot_marginal_rmst_curve

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:80, ]
x_cols <- grep("^x", names(dat), value = TRUE)[1:5]
X <- dat[, x_cols, drop = FALSE]
new.times <- seq(20, 120, by = 20)

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X,
  new.times = new.times,
  event.library = c("surv.coxph", "surv.glmnet"),
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE)
)

tau_grid <- seq(40, 120, by = 20)
plot_marginal_rmst_curve(
  fit = fit,
  data = dat,
  trt_col = "x4",
  times = new.times,
  tau_seq = tau_grid,
  inference = TRUE,
  B = 100,
  seed = 123
)




cleanEx()
nameEx("plot_patient_waterfall")
### * plot_patient_waterfall

flush(stderr()); flush(stdout())

### Name: plot_patient_waterfall
### Title: Waterfall Plot for an Individual Patient
### Aliases: plot_patient_waterfall

### ** Examples

if (requireNamespace("fastshap", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  shap_values <- explain_kernel(
    model = fit,
    X_explain = X[1:10, , drop = FALSE],
    X_background = X[11:40, , drop = FALSE],
    nsim = 5
  )

  plot_patient_waterfall(
    shap_values = shap_values,
    patient_index = 1,
    top_n = 5
  )
}



cleanEx()
nameEx("plot_predict")
### * plot_predict

flush(stderr()); flush(stdout())

### Name: plot_predict
### Title: Plot Predicted Survival Curves
### Aliases: plot_predict

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:120, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  eval_times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = eval_times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  preds <- predict(fit, newdata = X, new.times = eval_times, type = "event")

  plot_predict(
    preds = preds,
    eval_times = eval_times,
    patient_idx = c(1, 2)
  )
}



cleanEx()
nameEx("plot_rmst_vs_obs")
### * plot_rmst_vs_obs

flush(stderr()); flush(stdout())

### Name: plot_rmst_vs_obs
### Title: Plot Predicted RMST vs. Observed Survival Times
### Aliases: plot_rmst_vs_obs

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:80, ]
x_cols <- grep("^x", names(dat))[1:5]
X <- dat[, x_cols, drop = FALSE]
new.times <- seq(20, 120, by = 20)

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X,
  new.times = new.times,
  event.library = c("surv.coxph", "surv.glmnet"),
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE)
)

plot_rmst_vs_obs(
  fit = fit,
  data = dat,
  time_col = "duration",
  event_col = "event",
  times = new.times,
  tau = 350
)



cleanEx()
nameEx("plot_survival_heatmap")
### * plot_survival_heatmap

flush(stderr()); flush(stdout())

### Name: plot_survival_heatmap
### Title: Survival Probability Heatmap
### Aliases: plot_survival_heatmap

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  plot_survival_heatmap(
    object = fit,
    newdata = X[1:20, , drop = FALSE],
    times = times
  )
}



cleanEx()
nameEx("predict.SuperSurv")
### * predict.SuperSurv

flush(stderr()); flush(stdout())

### Name: predict.SuperSurv
### Title: Predict method for SuperSurv fits
### Aliases: predict.SuperSurv

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:80, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:10, , drop = FALSE]
  new.times <- seq(20, 120, by = 20)

  fit <- SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X,
    new.times = new.times,
    event.library = c("surv.coxph", "surv.ridge"),
    cens.library = c("surv.coxph"),
    control = list(saveFitLibrary = TRUE)
  )

  pred_event <- predict(
    object = fit,
    newdata = newX,
    new.times = new.times,
    type = "event"
  )

  dim(pred_event)
}



cleanEx()
nameEx("predict.surv.aorsf")
### * predict.surv.aorsf

flush(stderr()); flush(stdout())

### Name: predict.surv.aorsf
### Title: Prediction function for AORSF
### Aliases: predict.surv.aorsf

### ** Examples

if (requireNamespace("aorsf", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.aorsf(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    n_tree = 10,
    leaf_min_events = 2
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.bart")
### * predict.surv.bart

flush(stderr()); flush(stdout())

### Name: predict.surv.bart
### Title: Prediction function for BART
### Aliases: predict.surv.bart

### ** Examples

## Not run: 
##D if (.Platform$OS.type != "windows" &&
##D   requireNamespace("BART", quietly = TRUE)) {
##D   data("metabric", package = "SuperSurv")
##D   dat <- metabric[1:20, ]
##D   x_cols <- grep("^x", names(dat))[1:3]
##D   X <- dat[, x_cols, drop = FALSE]
##D   newX <- X[1:5, , drop = FALSE]
##D   times <- seq(50, 150, by = 50)
##D 
##D   fit <- surv.bart(
##D     time = dat$duration,
##D     event = dat$event,
##D     X = X,
##D     newdata = newX,
##D     new.times = times,
##D     ntree = 3,
##D     ndpost = 5,
##D     nskip = 5
##D   )
##D 
##D   pred <- fit[["pred"]]
##D   dim(pred)
##D }
## End(Not run)



cleanEx()
nameEx("predict.surv.coxboost")
### * predict.surv.coxboost

flush(stderr()); flush(stdout())

### Name: predict.surv.coxboost
### Title: Prediction function for CoxBoost wrapper
### Aliases: predict.surv.coxboost

### ** Examples

if (requireNamespace("CoxBoost", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.coxboost(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    stepno = 10,
    penalty = 50
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.coxph")
### * predict.surv.coxph

flush(stderr()); flush(stdout())

### Name: predict.surv.coxph
### Title: Prediction function for Cox regression wrapper
### Aliases: predict.surv.coxph

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.coxph(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
dim(pred)



cleanEx()
nameEx("predict.surv.gam")
### * predict.surv.gam

flush(stderr()); flush(stdout())

### Name: predict.surv.gam
### Title: Prediction function for GAM wrapper
### Aliases: predict.surv.gam

### ** Examples

if (requireNamespace("mgcv", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.gam(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    cts.num = 5
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.gbm")
### * predict.surv.gbm

flush(stderr()); flush(stdout())

### Name: predict.surv.gbm
### Title: Prediction function for GBM wrapper
### Aliases: predict.surv.gbm

### ** Examples

if (requireNamespace("gbm", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.gbm(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    n.trees = 20,
    interaction.depth = 1,
    shrinkage = 0.05,
    cv.folds = 0,
    n.minobsinnode = 3
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.glmnet")
### * predict.surv.glmnet

flush(stderr()); flush(stdout())

### Name: predict.surv.glmnet
### Title: Prediction function for GLMNET wrapper
### Aliases: predict.surv.glmnet

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.glmnet(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    alpha = 1,
    nfolds = 3
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.km")
### * predict.surv.km

flush(stderr()); flush(stdout())

### Name: predict.surv.km
### Title: Predict Method for Kaplan-Meier Wrapper
### Aliases: predict.surv.km

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.km(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
dim(pred)



cleanEx()
nameEx("predict.surv.parametric")
### * predict.surv.parametric

flush(stderr()); flush(stdout())

### Name: predict.surv.parametric
### Title: Prediction function for Universal Parametric Wrapper
### Aliases: predict.surv.parametric

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.parametric(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL,
  dist = "weibull"
)

pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
dim(pred)



cleanEx()
nameEx("predict.surv.ranger")
### * predict.surv.ranger

flush(stderr()); flush(stdout())

### Name: predict.surv.ranger
### Title: Prediction function for Ranger wrapper
### Aliases: predict.surv.ranger

### ** Examples

if (requireNamespace("ranger", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.ranger(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    num.trees = 10,
    min.node.size = 3
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.rfsrc")
### * predict.surv.rfsrc

flush(stderr()); flush(stdout())

### Name: predict.surv.rfsrc
### Title: Prediction function for survival random forest (RFSRC)
### Aliases: predict.surv.rfsrc

### ** Examples

if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.rfsrc(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    ntree = 10,
    nodesize = 3
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.ridge")
### * predict.surv.ridge

flush(stderr()); flush(stdout())

### Name: predict.surv.ridge
### Title: Prediction function for Ridge wrapper
### Aliases: predict.surv.ridge

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.ridge(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    nfolds = 3
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.rpart")
### * predict.surv.rpart

flush(stderr()); flush(stdout())

### Name: predict.surv.rpart
### Title: Prediction function for rpart wrapper
### Aliases: predict.surv.rpart

### ** Examples

if (requireNamespace("rpart", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.rpart(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    cp = 0.01,
    minsplit = 5,
    maxdepth = 3
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.svm")
### * predict.surv.svm

flush(stderr()); flush(stdout())

### Name: predict.surv.svm
### Title: Prediction function for SVM wrapper
### Aliases: predict.surv.svm

### ** Examples

if (requireNamespace("survivalsvm", quietly = TRUE) &&
  requireNamespace("quadprog", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:25, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.svm(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("predict.surv.xgboost")
### * predict.surv.xgboost

flush(stderr()); flush(stdout())

### Name: predict.surv.xgboost
### Title: Prediction function for XGBoost wrapper
### Aliases: predict.surv.xgboost

### ** Examples

if (requireNamespace("xgboost", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.xgboost(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    nrounds = 5,
    early_stopping_rounds = 2,
    max_depth = 1
  )

  pred <- predict(fit[["fit"]], newdata = newX, new.times = times)
  dim(pred)
}



cleanEx()
nameEx("screen.all")
### * screen.all

flush(stderr()); flush(stdout())

### Name: screen.all
### Title: Keep All Variables Screener
### Aliases: screen.all

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:20, ]
x_cols <- grep("^x", names(dat))[1:5]
X <- dat[, x_cols, drop = FALSE]

screen.all(X)



cleanEx()
nameEx("screen.elasticnet")
### * screen.elasticnet

flush(stderr()); flush(stdout())

### Name: screen.elasticnet
### Title: Elastic Net Screening Algorithm
### Aliases: screen.elasticnet

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]

  screen.elasticnet(
    time = dat$duration,
    event = dat$event,
    X = X,
    alpha = 0.5,
    minscreen = 2,
    nfolds = 3,
    nlambda = 20
  )
}



cleanEx()
nameEx("screen.glmnet")
### * screen.glmnet

flush(stderr()); flush(stdout())

### Name: screen.glmnet
### Title: GLMNET (Lasso) Screening
### Aliases: screen.glmnet

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]

  screen.glmnet(
    time = dat$duration,
    event = dat$event,
    X = X,
    alpha = 1,
    minscreen = 2,
    nfolds = 3,
    nlambda = 20
  )
}



cleanEx()
nameEx("screen.marg")
### * screen.marg

flush(stderr()); flush(stdout())

### Name: screen.marg
### Title: Marginal Cox Regression Screening
### Aliases: screen.marg

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:5]
X <- dat[, x_cols, drop = FALSE]

screen.marg(
  time = dat$duration,
  event = dat$event,
  X = X,
  minscreen = 2,
  min.p = 0.2
)



cleanEx()
nameEx("screen.rfsrc")
### * screen.rfsrc

flush(stderr()); flush(stdout())

### Name: screen.rfsrc
### Title: Random Survival Forest Screening Algorithm
### Aliases: screen.rfsrc

### ** Examples

if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]

  screen.rfsrc(
    time = dat$duration,
    event = dat$event,
    X = X,
    minscreen = 2,
    ntree = 10
  )
}



cleanEx()
nameEx("screen.var")
### * screen.var

flush(stderr()); flush(stdout())

### Name: screen.var
### Title: High Variance Screening Algorithm (Unsupervised)
### Aliases: screen.var

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:6]
X <- dat[, x_cols, drop = FALSE]

screen.var(
  time = dat$duration,
  event = dat$event,
  X = X,
  keep_fraction = 0.5,
  minscreen = 2
)



cleanEx()
nameEx("surv.aorsf")
### * surv.aorsf

flush(stderr()); flush(stdout())

### Name: surv.aorsf
### Title: Wrapper for AORSF (Oblique Random Survival Forest)
### Aliases: surv.aorsf

### ** Examples

if (requireNamespace("aorsf", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.aorsf(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    n_tree = 10,
    leaf_min_events = 2
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.bart")
### * surv.bart

flush(stderr()); flush(stdout())

### Name: surv.bart
### Title: Wrapper for BART (Bayesian Additive Regression Trees)
### Aliases: surv.bart

### ** Examples

## Not run: 
##D if (.Platform$OS.type != "windows" &&
##D   requireNamespace("BART", quietly = TRUE)) {
##D   data("metabric", package = "SuperSurv")
##D   dat <- metabric[1:20, ]
##D   x_cols <- grep("^x", names(dat))[1:3]
##D   X <- dat[, x_cols, drop = FALSE]
##D   newX <- X[1:5, , drop = FALSE]
##D   times <- seq(50, 150, by = 50)
##D 
##D   fit <- surv.bart(
##D     time = dat$duration,
##D     event = dat$event,
##D     X = X,
##D     newdata = newX,
##D     new.times = times,
##D     ntree = 3,
##D     ndpost = 5,
##D     nskip = 5
##D   )
##D 
##D   dim(fit[["pred"]])
##D }
## End(Not run)



cleanEx()
nameEx("surv.coxboost")
### * surv.coxboost

flush(stderr()); flush(stdout())

### Name: surv.coxboost
### Title: Wrapper function for Component-Wise Boosting (CoxBoost)
### Aliases: surv.coxboost

### ** Examples

if (requireNamespace("CoxBoost", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.coxboost(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    stepno = 10,
    penalty = 50
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.coxph")
### * surv.coxph

flush(stderr()); flush(stdout())

### Name: surv.coxph
### Title: Wrapper for standard Cox Proportional Hazards
### Aliases: surv.coxph

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.coxph(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.exponential")
### * surv.exponential

flush(stderr()); flush(stdout())

### Name: surv.exponential
### Title: Parametric Survival Prediction Wrapper (Exponential)
### Aliases: surv.exponential

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.exponential(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.gam")
### * surv.gam

flush(stderr()); flush(stdout())

### Name: surv.gam
### Title: Wrapper for Generalized Additive Cox Regression (GAM)
### Aliases: surv.gam

### ** Examples

if (requireNamespace("mgcv", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.gam(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    cts.num = 5
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.gbm")
### * surv.gbm

flush(stderr()); flush(stdout())

### Name: surv.gbm
### Title: Wrapper function for Gradient Boosting (GBM) prediction
###   algorithm
### Aliases: surv.gbm

### ** Examples

if (requireNamespace("gbm", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.gbm(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    n.trees = 20,
    interaction.depth = 1,
    shrinkage = 0.05,
    cv.folds = 0,
    n.minobsinnode = 3
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.glmnet")
### * surv.glmnet

flush(stderr()); flush(stdout())

### Name: surv.glmnet
### Title: Wrapper function for Penalized Cox Regression (GLMNET)
### Aliases: surv.glmnet

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.glmnet(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    alpha = 1,
    nfolds = 3
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.km")
### * surv.km

flush(stderr()); flush(stdout())

### Name: surv.km
### Title: Kaplan-Meier Prediction Algorithm
### Aliases: surv.km

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.km(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.loglogistic")
### * surv.loglogistic

flush(stderr()); flush(stdout())

### Name: surv.loglogistic
### Title: Parametric Survival Prediction Wrapper (Log-Logistic)
### Aliases: surv.loglogistic

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.loglogistic(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.lognormal")
### * surv.lognormal

flush(stderr()); flush(stdout())

### Name: surv.lognormal
### Title: Parametric Survival Prediction Wrapper (Log-Normal)
### Aliases: surv.lognormal

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.lognormal(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.parametric")
### * surv.parametric

flush(stderr()); flush(stdout())

### Name: surv.parametric
### Title: Universal Parametric Survival Wrapper
### Aliases: surv.parametric

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.parametric(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL,
  dist = "weibull"
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.ranger")
### * surv.ranger

flush(stderr()); flush(stdout())

### Name: surv.ranger
### Title: Wrapper function for Ranger Random Survival Forest
### Aliases: surv.ranger

### ** Examples

if (requireNamespace("ranger", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.ranger(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    num.trees = 10,
    min.node.size = 3
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.rfsrc")
### * surv.rfsrc

flush(stderr()); flush(stdout())

### Name: surv.rfsrc
### Title: Wrapper function for Random Survival Forests (RFSRC)
### Aliases: surv.rfsrc

### ** Examples

if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.rfsrc(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    ntree = 10,
    nodesize = 3
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.ridge")
### * surv.ridge

flush(stderr()); flush(stdout())

### Name: surv.ridge
### Title: Wrapper for Ridge Regression (Penalized Cox)
### Aliases: surv.ridge

### ** Examples

if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.ridge(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    nfolds = 3
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.rpart")
### * surv.rpart

flush(stderr()); flush(stdout())

### Name: surv.rpart
### Title: Wrapper for Survival Regression Trees (rpart)
### Aliases: surv.rpart

### ** Examples

if (requireNamespace("rpart", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.rpart(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    cp = 0.01,
    minsplit = 5,
    maxdepth = 3
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.svm")
### * surv.svm

flush(stderr()); flush(stdout())

### Name: surv.svm
### Title: Wrapper for Survival Support Vector Machine (survivalsvm)
### Aliases: surv.svm

### ** Examples

if (requireNamespace("survivalsvm", quietly = TRUE) &&
 requireNamespace("quadprog", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:25, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.svm(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times
  )

  dim(fit[["pred"]])
}



cleanEx()
nameEx("surv.weibull")
### * surv.weibull

flush(stderr()); flush(stdout())

### Name: surv.weibull
### Title: Parametric Survival Prediction Wrapper (Weibull)
### Aliases: surv.weibull

### ** Examples

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
newX <- X[1:5, , drop = FALSE]
times <- seq(50, 150, by = 50)

fit <- surv.weibull(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = newX,
  new.times = times,
  obsWeights = rep(1, nrow(dat)),
  id = NULL
)

dim(fit[["pred"]])



cleanEx()
nameEx("surv.xgboost")
### * surv.xgboost

flush(stderr()); flush(stdout())

### Name: surv.xgboost
### Title: Wrapper for XGBoost (Robust CV-Tuned + Safe Prediction)
### Aliases: surv.xgboost

### ** Examples

if (requireNamespace("xgboost", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:30, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  newX <- X[1:5, , drop = FALSE]
  times <- seq(50, 150, by = 50)

  fit <- surv.xgboost(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = newX,
    new.times = times,
    obsWeights = rep(1, nrow(dat)),
    id = NULL,
    nrounds = 5,
    early_stopping_rounds = 2,
    max_depth = 1
  )

  dim(fit[["pred"]])
}



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
