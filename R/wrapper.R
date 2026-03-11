#' Wrapper function for Random Survival Forests (RFSRC)
#'
#' Final Production Wrapper for RFSRC (Tunable & Robust).
#' Estimates a survival random forest using \code{\link[randomForestSRC]{rfsrc}}.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Currently ignored.
#' @param ntree Number of trees to grow (default: 1000).
#' @param nodesize Minimum number of deaths in terminal nodes (default: 15).
#' @param mtry Number of variables randomly selected as candidates for splitting a node.
#' @param ... Additional arguments passed to \code{\link[randomForestSRC]{rfsrc}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("randomForestSRC", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.rfsrc(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     ntree = 10,
#'     nodesize = 3
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.rfsrc <- function(time, event, X, newdata = NULL, new.times, obsWeights = NULL, id = NULL,
                       ntree = 1000, nodesize = 15, mtry = NULL, ...) {

  requireNamespace("randomForestSRC", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)
  Surv <- survival::Surv

  # 1. Handle Defaults
  if(is.null(obsWeights)) obsWeights <- rep(1, length(time))
  if(is.null(newdata)) newdata <- X

  # 2. Prepare Data (RFSRC strictly needs data.frames)
  X_df <- as.data.frame(X)
  newdata_df <- as.data.frame(newdata)
  data <- data.frame(time = time, event = event, X_df)

  # 3. Setup Arguments
  # FIX: Removed "survival::" prefix. Using as.formula protects the parser.
  args_list <- list(
    formula = as.formula("Surv(time, event) ~ ."),
    data = data,
    case.wt = obsWeights,
    ntree = ntree,
    nodesize = nodesize,
    importance = "none" # Faster fitting
  )

  # Safely append any additional arguments passed via ... (like from create_grid)
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    args_list <- c(args_list, extra_args)
  }

  # Only add mtry if explicitly provided, otherwise let rfsrc calculate default
  if (!is.null(mtry)) args_list$mtry <- mtry

  # 4. Fit Model
  fit.rfsrc <- do.call(randomForestSRC::rfsrc, args_list)

  # 5. Predict on newdata (Training Grid)
  survs <- predict(fit.rfsrc, newdata = newdata_df, importance = "none")$survival
  train_times <- fit.rfsrc$time.interest

  # 6. Interpolate to new.times
  # Uses constant interpolation for step-function survival curves
  pred <- t(apply(survs, 1, function(y) {
    stats::approx(train_times, y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  if (is.null(dim(pred))) {
    pred <- matrix(pred, nrow = nrow(newdata_df), ncol = length(new.times))
  }

  # 7. Safety: Monotonicity and Clamping
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    # Ensure survival curve never goes up
    pred <- t(apply(pred, 1, cummin))
  }

  fit <- list(object = fit.rfsrc)
  class(fit) <- c("surv.rfsrc")
  list(pred = pred, fit = fit)
}



#' Prediction function for survival random forest (RFSRC)
#'
#' Obtains predicted survivals from a fitted \code{surv.rfsrc} object.
#'
#' @param object Fitted \code{surv.rfsrc} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("randomForestSRC", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.rfsrc(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     ntree = 10,
#'     nodesize = 3
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
# 1. Change the argument here to 'newdata'
predict.surv.rfsrc <- function(object, newdata, new.times, ...) {


  fit <- object$object

  # 3. Predict on native time grid (ensure you use newdata here!)
  survs <- predict(fit, newdata = as.data.frame(newdata), importance = "none")$survival
  train_times <- fit$time.interest

  # 4. Vectorized interpolation to new.times
  pred <- t(apply(survs, 1, function(y) {
    stats::approx(train_times, y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  # 5. Safety Check: Dimensions
  if (is.null(dim(pred))) {
    pred <- matrix(pred, nrow = nrow(newdata), ncol = length(new.times))
  }

  # 6. Safety Check: Clamp and Monotonicity
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}




############################################################
## XGBoost (Cox) — survSuperLearner wrapper + predict method
############################################################

#' Wrapper for XGBoost (Robust CV-Tuned + Safe Prediction)
#'
#' Estimates a Cox proportional hazards model via XGBoost.
#' Incorporates safe Breslow hazard calculation and matrix alignment to prevent C++ crashes.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param nrounds Max number of boosting iterations (default: 1000).
#' @param early_stopping_rounds Rounds with no improvement to trigger early stopping (default: 10).
#' @param eta Learning rate (default: 0.05).
#' @param max_depth Maximum tree depth (default: 2).
#' @param min_child_weight Minimum sum of instance weight in a child (default: 5).
#' @param lambda L2 regularization term on weights (default: 10).
#' @param subsample Subsample ratio of the training instances (default: 0.7).
#' @param ... Additional arguments passed to \code{\link[xgboost]{xgb.train}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("xgboost", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.xgboost(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     nrounds = 5,
#'     early_stopping_rounds = 2,
#'     max_depth = 1
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.xgboost <- function(time, event, X, newdata = NULL,  new.times, obsWeights, id,
                         nrounds = 1000, early_stopping_rounds = 10,
                         eta = 0.05, max_depth = 2, min_child_weight = 5,
                         lambda = 10, subsample = 0.7, ...) {

  requireNamespace("xgboost", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  if(is.null(newdata)) newdata <- X

  # 1. Prepare Matrices
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  # CRITICAL FIX: Matrix Alignment for CV
  # If newdata is missing a factor level present in X, we must pad it with zeros
  missing_cols <- setdiff(colnames(X_mat), colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  # Ensure columns are in the exact same order
  newdata_mat <- newdata_mat[, colnames(X_mat), drop = FALSE]

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 2. Format for XGBoost Cox
  y_xgb <- ifelse(event == 1, time, -time)
  dtrain <- xgboost::xgb.DMatrix(data = X_mat, label = y_xgb, weight = obsWeights)

  # 3. Fit XGBoost
  fit <- xgboost::xgb.train(
    params = list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = eta,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      lambda = lambda,
      subsample = subsample
    ),
    data = dtrain,
    nrounds = nrounds,
    watchlist = list(train = dtrain),
    early_stopping_rounds = early_stopping_rounds,
    verbose = 0,
    ...
  )

  # 4. Center Training Scores (To stabilize hazard)
  lp_train <- predict(fit, newdata = X_mat)
  tr_mean <- mean(lp_train)
  lp_train_centered <- lp_train - tr_mean

  # 5. Calculate Baseline Hazard (using your custom safe function)
  bh <- safe_breslow_step(
    time = time, event = event,
    risk_score = lp_train_centered, new.times = new.times
  )

  # 6. Predict on New Data
  lp_new <- predict(fit, newdata = newdata_mat)
  lp_new_centered <- lp_new - tr_mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 7. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1


  fit_obj <- list(
    object  = fit,
    basehaz = bh,
    times   = new.times,
    stats   = list(mean = tr_mean, features = colnames(X_mat))
  )
  class(fit_obj) <- c("surv.xgboost")


  list(pred = pred, fit = fit_obj)
}





#' Prediction function for XGBoost wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.xgboost} object.
#'
#' @param object Fitted \code{surv.xgboost} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("xgboost", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.xgboost(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     nrounds = 5,
#'     early_stopping_rounds = 2,
#'     max_depth = 1
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.xgboost <- function(object, newdata, new.times, ...) {


  # 1. Faster/Safer Interpolation using stepfun
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Format Matrix and Align Columns
  # Must use the exact features the XGBoost model expects
  expected_cols <- object$stats$features
  if (is.null(expected_cols)) {
    expected_cols <- object$object$feature_names
  }
  if (is.null(expected_cols)) stop("Cannot determine expected feature names for xgboost model.")


  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(expected_cols, colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, expected_cols, drop = FALSE]

  # 3. Predict Linear Predictor
  lp_new <- predict(object$object, newdata = newdata_mat)

  # 4. Center and Clamp
  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 5. Survival Formula
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 6. Safety Clamping
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}







############################################################
## Survival SVM — Wrapper
## Pattern: Utility Score -> Negate -> Cox Offset -> Survival
############################################################

#' Wrapper for Survival Support Vector Machine (survivalsvm)
#'
#' Final Production Wrapper for SVM (Tunable & Robust).
#' Estimates a survival SVM and calibrates the raw utility scores into
#' survival probabilities using a univariate Cox proportional hazards model.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param gamma.mu Regularization parameter for the SVM (default: 0.1).
#' @param type Type of SVM implementation (default: "vanbelle2").
#' @param kernel Kernel type for the SVM (default: "lin_kernel").
#' @param opt.meth Optimization method (default: "quadprog").
#' @param ... Additional arguments passed to \code{\link[survivalsvm]{survivalsvm}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("survivalsvm", quietly = TRUE) &&
#'  requireNamespace("quadprog", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:25, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.svm(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.svm <- function(time, event, X, newdata, new.times, obsWeights, id,
                     gamma.mu = 0.1, type = "vanbelle2",
                     kernel = "lin_kernel", opt.meth = "quadprog", ...) {

  requireNamespace("survivalsvm", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Prepare Data
  Y <- survival::Surv(time, event)
  df_train <- data.frame(Y = Y, X)

  # 2. Fit Survival SVM
  makediff3 <- utils::getFromNamespace("makediff3", "survivalsvm")
  fit <- survivalsvm::survivalsvm(
    formula = Y ~ .,
    data = df_train,
    type = type,
    diff.meth = "makediff3",
    gamma.mu = gamma.mu,
    kernel = kernel,
    opt.meth = opt.meth,
    ...
  )

  # 3. Extract Raw Utility Scores and Negate them to represent Risk
  p_tr <- predict(fit, newdata = X)
  lp_raw <- -1 * as.numeric(p_tr$predicted)

  # 4. Fit Calibrator Cox Model
  calib_data <- data.frame(time = time, event = event, score = lp_raw)
  calib_fit <- survival::coxph(
    survival::Surv(time, event) ~ score,
    data = calib_data,
    weights = obsWeights
  )

  # 5. Extract Baseline Hazard from Calibrator
  bh_df <- survival::basehaz(calib_fit, centered = FALSE)
  bh <- stats::approx(
    bh_df$time, bh_df$hazard, xout = new.times,
    method = "constant", f = 0, rule = 2, ties = mean
  )$y
  bh <- cummax(replace(bh, is.na(bh), 0))

  # 6. Predict and Calibrate on New Data
  p_te <- predict(fit, newdata = newdata)
  lp_new_raw <- -1 * as.numeric(p_te$predicted)

  # Use R's native predict.coxph to get perfectly centered linear predictors
  lp_new_calibrated <- predict(
    calib_fit,
    newdata = data.frame(score = lp_new_raw),
    type = "lp"
  )

  lp_new_calibrated <- pmin(lp_new_calibrated, 700)

  # 7. Convert to Survival Probabilities
  pred <- outer(lp_new_calibrated, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(
    object = fit,
    calibrator = calib_fit,
    basehaz = bh,
    times = new.times
  )
  class(fit_obj) <- c("surv.svm")

  list(pred = pred, fit = fit_obj)
}

#' Prediction function for SVM wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.svm} object.
#' Applies the learned Cox calibration to the new raw utility scores.
#'
#' @param object Fitted \code{surv.svm} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("survivalsvm", quietly = TRUE) &&
#'   requireNamespace("quadprog", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:25, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.svm(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.svm <- function(object, newdata, new.times, ...) {

  # 1. Faster/Safer Interpolation using stepfun
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Get Raw Predictions from SVM
  p_te <- predict(object$object, newdata = newdata)
  lp_new_raw <- -1 * as.numeric(p_te$predicted)

  # 3. Apply Calibration using the saved Cox model
  lp_new_calibrated <- predict(
    object$calibrator,
    newdata = data.frame(score = lp_new_raw),
    type = "lp"
  )

  # 4. Convert to Survival Probabilities S(t)
  pred <- outer(lp_new_calibrated, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 5. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  return(pred)
}









#' Wrapper for Survival Regression Trees (rpart)
#'
#' Final Production Wrapper for single decision trees.
#' Uses \code{\link[rpart]{rpart}} with method="exp" and calculates
#' survival probabilities using the Breslow estimator.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param cp Complexity parameter (default: 0.01).
#' @param minsplit Minimum number of observations to attempt a split (default: 20).
#' @param maxdepth Maximum depth of any node of the final tree (default: 30).
#' @param ... Additional arguments passed to \code{\link[rpart]{rpart.control}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.rpart(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     cp = 0.01,
#'     minsplit = 5,
#'     maxdepth = 3
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.rpart <- function(time, event, X, newdata, new.times, obsWeights, id,
                       cp = 0.01, minsplit = 20, maxdepth = 30, ...) {

  requireNamespace("rpart", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Safety Clamp for Time
  time <- pmax(time, 1e-5)

  # 2. Prepare Data
  X_df <- as.data.frame(X)
  newdata_df <- as.data.frame(newdata)
  dat <- data.frame(time = time, event = event, X_df)

  # 3. Fit the Survival Tree
  fit <- rpart::rpart(
    survival::Surv(time, event) ~ .,
    data = dat,
    weights = obsWeights,
    method = "exp", # Survival method for rpart
    control = rpart::rpart.control(
      cp = cp,
      minsplit = minsplit,
      maxdepth = maxdepth,
      ...
    )
  )

  # 4. Predict Expected Event Rates and Convert to Linear Predictor
  # If a leaf has 0 events, the rate is 0. log(0) = -Inf, which crashes Breslow.
  # We use pmax to safely clamp the minimum rate before taking the log.
  rate_train <- predict(fit, newdata = X_df)
  lp_train <- log(pmax(rate_train, 1e-10))

  # 5. Center Training Scores
  tr_mean <- mean(lp_train)
  lp_train_centered <- lp_train - tr_mean

  # 6. Baseline Hazard Calculation
  bh <- safe_breslow_step(
    time = time,
    event = event,
    risk_score = lp_train_centered,
    new.times = new.times
  )

  # 7. Predict on newdata
  rate_new <- predict(fit, newdata = newdata_df)
  lp_new <- log(pmax(rate_new, 1e-10))

  lp_new_centered <- lp_new - tr_mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 8. Calculate Survival S(t)
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata_df), ncol = length(new.times))

  # 9. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  fit_obj <- list(
    object = fit,
    basehaz = bh,
    times = new.times,
    stats = list(mean = tr_mean)
  )
  class(fit_obj) <- c("surv.rpart")

  list(pred = pred, fit = fit_obj)
}


#' Prediction function for rpart wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.rpart} object.
#'
#' @param object Fitted \code{surv.rpart} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.rpart(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     cp = 0.01,
#'     minsplit = 5,
#'     maxdepth = 3
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.rpart <- function(object, newdata, new.times, ...) {

  # 1. Reconstruct Baseline Hazard
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Predict Event Rate and Convert to Centered LP
  rate_new <- predict(object$object, newdata = as.data.frame(newdata))
  lp_new <- log(pmax(rate_new, 1e-10))

  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 3. Calculate Survival S(t)
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 4. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}











#' Wrapper for Ridge Regression (Penalized Cox)
#'
#' Final Production Wrapper for Ridge Regression (Tunable & Robust).
#' Estimates a penalized Cox model using a pure Ridge penalty (alpha = 0).
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param nfolds Number of folds for internal cross-validation to select lambda. Default is 10.
#' @param ... Additional arguments passed to \code{\link[glmnet]{cv.glmnet}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.ridge(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     nfolds = 3
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.ridge <- function(time, event, X, newdata, new.times, obsWeights = NULL, id = NULL, nfolds = 10, ...) {

  requireNamespace("glmnet", quietly = TRUE)

  # We clamp the minimum time to a tiny positive number.
  time <- pmax(time, 1e-5)

  if (is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Prepare Matrices & Align Columns
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(colnames(X_mat), colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, colnames(X_mat), drop = FALSE]

  # 2. Fit Ridge (Hardcoded alpha = 0 for pure Ridge)
  fit <- glmnet::cv.glmnet(
    x = X_mat,
    y = survival::Surv(time, event),
    family = "cox",
    alpha = 0,               # Enforce Ridge penalty
    weights = obsWeights,
    nfolds = nfolds,
    ...
  )

  # 3. Predict Risk Scores & Center them
  lp_train <- as.numeric(predict(fit, newx = X_mat, s = "lambda.min", type = "link"))
  tr_mean <- mean(lp_train)
  lp_train_centered <- lp_train - tr_mean

  # 4. Baseline Hazard (Using our robust helper)
  bh <- safe_breslow_step(
    time = time,
    event = event,
    risk_score = lp_train_centered,
    new.times = new.times
  )

  # 5. Predict on newdata
  lp_new <- as.numeric(predict(fit, newx = newdata_mat, s = "lambda.min", type = "link"))
  lp_new_centered <- lp_new - tr_mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(
    object = fit,
    basehaz = bh,
    times = new.times,
    stats = list(mean = tr_mean, features = colnames(X_mat))
  )
  class(fit_obj) <- c("surv.ridge")

  list(pred = pred, fit = fit_obj)
}


#' Prediction function for Ridge wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.ridge} object.
#'
#' @param object Fitted \code{surv.ridge} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.ridge(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     nfolds = 3
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.ridge <- function(object, newdata, new.times, ...) {

  # 1. Faster/Safer Interpolation
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Format Matrix and Align Columns
  expected_cols <- object$stats$features
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(expected_cols, colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, expected_cols, drop = FALSE]

  # 3. Predict Linear Predictor
  lp_new <- as.numeric(predict(object$object, newx = newdata_mat, s = "lambda.min", type = "link"))

  # 4. Center and Clamp
  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 5. Survival Formula
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}






############################################################
## Ranger (Random Forest) — wrapper + predict method
## Pattern: Native Probability Output + Interpolation
############################################################

#' Wrapper function for Ranger Random Survival Forest
#'
#' Final Production Wrapper for Ranger (Tunable & Fast).
#' Uses the \code{\link[ranger]{ranger}} C++ implementation to estimate survival curves.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param num.trees Number of trees (default: 500).
#' @param mtry Number of variables to split at each node. Defaults to \code{sqrt(p)}.
#' @param min.node.size Minimum node size (default: 15 for survival).
#' @param ... Additional arguments passed to \code{\link[ranger]{ranger}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.ranger(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     num.trees = 10,
#'     min.node.size = 3
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.ranger <- function(time, event, X, newdata, new.times, obsWeights, id,
                        num.trees = 500, mtry = NULL, min.node.size = NULL, ...) {

  requireNamespace("ranger", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Prepare Data
  dat <- data.frame(time = time, status = event, as.data.frame(X))

  # 2. Setup Arguments
  args_list <- list(
    formula = survival::Surv(time, status) ~ .,
    data = dat,
    num.trees = num.trees,
    case.weights = obsWeights,
    probability = FALSE, # Survival uses probability=FALSE in ranger
    verbose = FALSE,
    ...
  )

  if (!is.null(mtry)) args_list$mtry <- mtry
  if (!is.null(min.node.size)) args_list$min.node.size <- min.node.size

  # 3. Fit Model
  fit <- do.call(ranger::ranger, args_list)

  # 4. Predict on newdata (Training Grid)
  p_obj <- predict(fit, data = as.data.frame(newdata))
  surv_probs <- p_obj$survival
  train_times <- p_obj$unique.death.times

  # 5. Vectorized Interpolation to new.times
  pred <- t(apply(surv_probs, 1, function(y) {
    stats::approx(x = train_times, y = y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  # 6. Safety Clamp and Monotonicity
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  fit_obj <- list(object = fit, times = new.times)
  class(fit_obj) <- c("surv.ranger")

  list(pred = pred, fit = fit_obj)
}




#' Prediction function for Ranger wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.ranger} object.
#'
#' @param object Fitted \code{surv.ranger} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.ranger(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     num.trees = 10,
#'     min.node.size = 3
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.ranger <- function(object, newdata, new.times, ...) {

  # 1. Predict using saved model
  p_obj <- predict(object$object, data = as.data.frame(newdata))
  surv_probs  <- p_obj$survival
  train_times <- p_obj$unique.death.times

  # 2. Vectorized Interpolation
  pred <- t(apply(surv_probs, 1, function(y) {
    stats::approx(x = train_times, y = y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  # 3. Dimensions Safety Check
  if (is.null(dim(pred))) {
    pred <- matrix(pred, nrow = nrow(newdata), ncol = length(new.times))
  }

  # 4. Safety Clamp and Monotonicity
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}







#' Universal Parametric Survival Wrapper
#'
#' Final Production Wrapper for AFT Models (Weibull, Exponential, LogNormal, LogLogistic).
#' Replaces individual wrappers with one robust, vectorized function.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param dist Distribution for the AFT model (default: "weibull").
#' @param ... Additional arguments passed to \code{\link[survival]{survreg}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.parametric(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL,
#'   dist = "weibull"
#' )
#'
#' dim(fit$pred)
#' @export
surv.parametric <- function(time, event, X, newdata, new.times, obsWeights, id, dist = "weibull", ...) {
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Handle Time <= 0 issue for AFT models (log(0) = -Inf)
  time <- pmax(time, 1e-5)

  # 2. Fit Model
  # CRITICAL FIX: Use data.frame() to prevent cbind from destroying factor variables
  dat <- data.frame(time = time, event = event, as.data.frame(X))

  fit <- survival::survreg(
    survival::Surv(time, event) ~ .,
    data = dat,
    weights = obsWeights,
    dist = dist,
    ...
  )

  # 3. Fast Prediction (Closed Form)
  # predict.survreg uses type = "linear" to return the linear predictor
  lp <- predict(fit, newdata = as.data.frame(newdata), type = "linear")
  scale <- fit$scale

  # 4. Vectorized Survival Calculation
  if (dist %in% c("exponential", "weibull")) {
    pred <- outer(lp, new.times, function(eta, t) {
      exp( - (t / exp(eta))^(1/scale) )
    })

  } else if (dist == "lognormal") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 - stats::pnorm((log(t) - eta) / scale)
    })

  } else if (dist == "loglogistic") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 / (1 + (t / exp(eta))^(1/scale))
    })

  } else {
    stop("Distribution not supported in fast wrapper")
  }

  # 5. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(object = fit, dist = dist)
  class(fit_obj) <- c("surv.parametric")

  list(pred = pred, fit = fit_obj)
}





#' Parametric Survival Prediction Wrapper (Exponential)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.exponential(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' dim(fit$pred)
#' @export
surv.exponential <- function(time, event, X, newdata, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newdata = newdata,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "exponential", ...)
}



#' Parametric Survival Prediction Wrapper (Log-Logistic)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.loglogistic(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' dim(fit$pred)
#' @export
surv.loglogistic <- function(time, event, X, newdata, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newdata = newdata,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "loglogistic", ...)
}



#' Parametric Survival Prediction Wrapper (Log-Normal)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.lognormal(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' dim(fit$pred)
#' @export
surv.lognormal <- function(time, event, X, newdata, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newdata = newdata,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "lognormal", ...)
}

#' Parametric Survival Prediction Wrapper (Weibull)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.weibull(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' dim(fit$pred)
#' @export
surv.weibull <- function(time, event, X, newdata, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newdata = newdata,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "weibull", ...)
}






#' Prediction function for Universal Parametric Wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.parametric} object
#' using exact closed-form equations.
#'
#' @param object Fitted \code{surv.parametric} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.parametric(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL,
#'   dist = "weibull"
#' )
#'
#' pred <- predict(fit$fit, newdata = newX, new.times = times)
#' dim(pred)
#' @export
predict.surv.parametric <- function(object, newdata, new.times, ...) {

  fit <- object$object
  dist <- object$dist
  scale <- fit$scale

  # Extract linear predictor
  lp <- predict(fit, newdata = as.data.frame(newdata), type = "linear")

  # Vectorized Survival Calculation
  if (dist %in% c("exponential", "weibull")) {
    pred <- outer(lp, new.times, function(eta, t) {
      exp( - (t / exp(eta))^(1/scale) )
    })
  } else if (dist == "lognormal") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 - stats::pnorm((log(t) - eta) / scale)
    })
  } else if (dist == "loglogistic") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 / (1 + (t / exp(eta))^(1/scale))
    })
  }

  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  return(pred)
}







#' Kaplan-Meier Prediction Algorithm
#'
#' This prediction algorithm ignores all covariates and computes the marginal
#' Kaplan-Meier survival estimator using the \code{\link[survival]{survfit}} function.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param X Training covariate data.frame (Ignored by KM).
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Numeric vector of times at which to predict survival.
#' @param obsWeights Numeric vector of observation weights.
#' @param id Optional vector indicating subject/cluster identities.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: A list containing the fitted \code{\link[survival]{survfit}} object.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at \code{new.times}.
#' }
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.km(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' dim(fit$pred)
#' @export
surv.km <- function(time, event, X, newdata, new.times, obsWeights, id, ...) {

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # Fit the weighted KM curve
  fit.km <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    weights = obsWeights
  )

  # Efficient prediction using stepfun
  sfun <- stats::stepfun(fit.km$time, c(1, fit.km$surv), right = FALSE)
  surv_probs <- sfun(new.times)

  # Repeat for every patient (KM gives the exact same curve for everyone)
  pred <- matrix(surv_probs,
                 nrow = nrow(newdata),
                 ncol = length(new.times),
                 byrow = TRUE)

  fit <- list(object = fit.km)
  class(fit) <- c("surv.km")

  return(list(pred = pred, fit = fit))
}



#' Predict Method for Kaplan-Meier Wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.km} object.
#'
#' @param object A fitted object of class \code{surv.km}.
#' @param newdata New covariate data.frame for which to obtain predictions (Ignored).
#' @param new.times Numeric vector of times at which to predict survival.
#' @param ... Additional ignored arguments.
#'
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation times.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.km(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' pred <- predict(fit$fit, newdata = newX, new.times = times)
#' dim(pred)
#' @export
predict.surv.km <- function(object, newdata, new.times, ...) {

  # Extract the original fitted KM object
  fit.km <- object$object

  # Reconstruct the step function
  sfun <- stats::stepfun(fit.km$time, c(1, fit.km$surv), right = FALSE)
  surv_probs <- sfun(new.times)

  # Replicate across all new patients
  pred_matrix <- matrix(surv_probs,
                        nrow = nrow(newdata),
                        ncol = length(new.times),
                        byrow = TRUE)

  return(pred_matrix)
}








#' Wrapper function for Penalized Cox Regression (GLMNET)
#'
#' Final Production Wrapper for GLMNET (Tunable & Robust).
#' Estimates a penalized Cox model (Lasso, Ridge, or Elastic Net) with automatic lambda selection.
#' Uses the Breslow estimator with a step-function approach for the baseline hazard.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param alpha The elasticnet mixing parameter (0 = Ridge, 1 = Lasso). Default is 1.
#' @param nfolds Number of folds for internal cross-validation to select lambda. Default is 10.
#' @param ... Additional arguments passed to \code{\link[glmnet]{cv.glmnet}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }

#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.glmnet(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     alpha = 1,
#'     nfolds = 3
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.glmnet <- function(time, event, X, newdata, new.times, obsWeights, id,
                        alpha = 1, nfolds = 10, ...) {

  requireNamespace("glmnet", quietly = TRUE)

  # We clamp the minimum time to a tiny positive number.
  time <- pmax(time, 1e-5)

  # 1. Convert factors to dummy variables (glmnet requires matrices)
  X_mat <- stats::model.matrix(~ . - 1, data = X)
  newdata_mat <- stats::model.matrix(~ . - 1, data = newdata)

  missing_cols <- setdiff(colnames(X_mat), colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, colnames(X_mat), drop = FALSE]

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 2. Fit the Penalized Model
  # We use standard argument passing (cleaner than do.call here)
  fit <- glmnet::cv.glmnet(
    x = X_mat,
    y = survival::Surv(time, event),
    family = "cox",
    weights = obsWeights,
    alpha = alpha,
    nfolds = nfolds,
    ...
  )

  beta_min <- as.matrix(stats::coef(fit, s = "lambda.min"))
  s_use <- "lambda.min"

  if (sum(beta_min != 0) == 0) {
    beta_path <- as.matrix(fit$glmnet.fit$beta)
    nnz_path <- colSums(beta_path != 0)

    idx_nonzero <- which(nnz_path > 0)
    if (length(idx_nonzero) > 0) {
      s_use <- fit$glmnet.fit$lambda[idx_nonzero[1]]
    }
  }

  # 3. Calibrate Baseline Hazard (Breslow Estimator)
  lp_train <- as.numeric(predict(fit, newx = X_mat, s = s_use, type = "link"))
  cox_off <- survival::coxph(
    survival::Surv(time, event) ~ offset(lp_train),
    weights = obsWeights,
    ties = "breslow"
  )

  bh_df <- survival::basehaz(cox_off, centered = FALSE)
  bh <- stats::approx(
    bh_df$time, bh_df$hazard, xout = new.times,
    method = "constant", f = 0, rule = 2, ties = mean
  )$y
  bh <- cummax(replace(bh, is.na(bh), 0))

  # 4. Predict on newdata
  lp_new <- as.numeric(predict(fit, newx = newdata_mat, s = s_use, type = "link"))
  pred <- outer(lp_new, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 5. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(
    object = fit,
    basehaz = bh,
    times = new.times,
    stats = list(
      features = colnames(X_mat),
      s_use = s_use
    )
  )
  class(fit_obj) <- c("surv.glmnet")

  list(pred = pred, fit = fit_obj)
}




#' Prediction function for GLMNET wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.glmnet} object.
#'
#' @param object Fitted \code{surv.glmnet} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.glmnet(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     alpha = 1,
#'     nfolds = 3
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.glmnet <- function(object, newdata, new.times, ...) {

  # 1. Align Baseline Hazard (Step Function)
  if (identical(all.equal(new.times, object$times), TRUE)) {
    bh <- object$basehaz
  } else {
    bh <- stats::approx(
      object$times, object$basehaz, xout = new.times,
      method = "constant", f = 0, rule = 2, ties = mean
    )$y
  }
  bh <- cummax(replace(bh, is.na(bh), 0))

  # 2. Convert newdata to matrix safely
  # model.matrix removes NAs, so we ensure the dataframe is intact

  expected_cols <- object$stats$features
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(expected_cols, colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, expected_cols, drop = FALSE]


  # 3. Predict Linear Predictor
  s_use <- object$stats$s_use
  if (is.null(s_use)) s_use <- "lambda.min"

  lp_new <- as.numeric(predict(
    object$object,
    newx = newdata_mat,
    s = s_use,
    type = "link"
  ))


  # 4. Convert LP to Survival Probability Matrix
  pred <- outer(lp_new, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 5. Safety Clamp
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}






#' Wrapper function for Gradient Boosting (GBM) prediction algorithm
#'
#' Final Production Wrapper for GBM (Tunable & Robust).
#' Estimates a Cox proportional hazards model via gradient boosting.
#' Uses the Breslow estimator with a step-function approach for the baseline hazard.
#' Includes internal safeguards against C++ crashes and small cross-validation folds.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param n.trees Integer specifying the total number of trees to fit (default: 1000).
#' @param interaction.depth Maximum depth of variable interactions (default: 2).
#' @param shrinkage A shrinkage parameter applied to each tree (default: 0.01).
#' @param cv.folds Number of cross-validation folds to perform internally for optimal tree selection (default: 5).
#' @param n.minobsinnode Minimum number of observations in the trees terminal nodes (default: 10).
#' @param ... Additional arguments passed to \code{\link[gbm]{gbm}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("gbm", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.gbm(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     n.trees = 20,
#'     interaction.depth = 1,
#'     shrinkage = 0.05,
#'     cv.folds = 0,
#'     n.minobsinnode = 3
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.gbm <- function(time, event, X, newdata, new.times, obsWeights, id,
                     n.trees = 1000, interaction.depth = 2, shrinkage = 0.01,
                     cv.folds = 5, n.minobsinnode = 10, ...) {

  requireNamespace("gbm", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  # 1. CRITICAL FIX: Prevent C++ Crashes (Data Types)
  if(is.matrix(X)) X <- as.data.frame(X)
  if(is.matrix(newdata)) newdata <- as.data.frame(newdata)

  # Force characters to factors
  X[] <- lapply(X, function(x) if(is.character(x)) as.factor(x) else x)
  newdata[] <- lapply(newdata, function(x) if(is.character(x)) as.factor(x) else x)

  # Align factor levels in newdata to match training
  for (nm in intersect(names(X), names(newdata))) {
    if (is.factor(X[[nm]])) {
      newdata[[nm]] <- factor(newdata[[nm]], levels = levels(X[[nm]]))
    }
  }


  # 2. CRITICAL FIX: Safety for Small Folds
  # If a CV fold is tiny, n.minobsinnode=10 will crash GBM.
  n_obs <- nrow(X)
  safe_min_node <- max(3, floor(n_obs * 0.05)) # Ensure at least 3 obs per node
  n.minobsinnode <- min(n.minobsinnode, safe_min_node)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))
  dat <- data.frame(time = time, status = event, X)

  # 3. Fit GBM Model
  # We use explicit arguments so create_surv_grid can tune them perfectly
  fit <- gbm::gbm(
    formula = survival::Surv(time, status) ~ .,
    data = dat,
    distribution = "coxph",
    n.trees = n.trees,
    interaction.depth = interaction.depth,
    shrinkage = shrinkage,
    cv.folds = cv.folds,
    n.minobsinnode = n.minobsinnode,
    weights = obsWeights,
    verbose = FALSE,
    keep.data = FALSE,
    ...
  )

  # 4. Get Best Iteration
  method_perf <- if(cv.folds > 0) "cv" else "OOB"
  best.iter <- gbm::gbm.perf(fit, method = method_perf, plot.it = FALSE)

  # 5. Calibration (Breslow Estimator)
  # Get Risk Scores (Linear Predictor) on Training Data
  lp_train <- predict(fit, newdata = X, n.trees = best.iter, type = "link")

  cox_df <- data.frame(time = time, event = event, lp = lp_train)
  cox_off <- survival::coxph(
    survival::Surv(time, event) ~ offset(lp),
    data = cox_df,
    weights = obsWeights,
    ties = "breslow"
  )

  bh_df <- survival::basehaz(cox_off, centered = FALSE)

  # Interpolate to new.times (Step Function)
  bh <- stats::approx(
    bh_df$time, bh_df$hazard, xout = new.times,
    method = "constant", f = 0, rule = 2, ties = mean
  )$y
  bh[is.na(bh)] <- 0
  bh <- cummax(bh)

  # 6. Predict on New Data
  lp_new <- predict(fit, newdata = newdata, n.trees = best.iter, type = "link")

  # S(t) = exp(-H0(t) * exp(lp))
  pred <- outer(lp_new, bh, function(lp, h) exp(-exp(lp) * h))

  # 7. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(object = fit, basehaz = bh, times = new.times, best.iter = best.iter)
  class(fit_obj) <- c("surv.gbm")

  list(pred = pred, fit = fit_obj)
}


#' Prediction function for GBM wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.gbm} object.
#' Uses a step function to align the Breslow baseline hazard with the requested times.
#'
#' @param object Fitted \code{surv.gbm} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("gbm", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.gbm(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     n.trees = 20,
#'     interaction.depth = 1,
#'     shrinkage = 0.05,
#'     cv.folds = 0,
#'     n.minobsinnode = 3
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.gbm <- function(object, newdata, new.times, ...) {

  # 1. Safer Interpolation using stepfun
  # Prepend 0 to handle times before the first event safely
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Linear Predictor
  lp_new <- predict(
    object$object,
    newdata = as.data.frame(newdata),
    n.trees = object$best.iter,
    type = "link"
  )

  # 3. Survival Formula S(t) = exp(-H0(t) * exp(lp))
  pred <- outer(lp_new, bh, function(lp, h) exp(-exp(lp) * h))

  # 4. Safety Clamping
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}








#' Wrapper for Generalized Additive Cox Regression (GAM)
#'
#' Final Production Wrapper for GAM (Tunable & Robust).
#' Uses \code{\link[mgcv]{gam}} to fit an additive combination of smooth
#' and linear functions.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights (Note: Ignored, as mgcv uses weights for the event indicator).
#' @param id Optional cluster/individual ID indicator.
#' @param cts.num Cutoff of unique values at which a numeric covariate receives a smooth term (s).
#' @param ... Additional arguments passed to \code{\link[mgcv]{gam}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.gam(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     cts.num = 5
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.gam <- function(time, event, X, newdata, new.times, obsWeights, id, cts.num = 5, ...) {

  requireNamespace("mgcv", quietly = TRUE)

  # 1. Formula Construction (Safe Spline Assignment)
  X_df <- as.data.frame(X)
  newdata_df <- as.data.frame(newdata)

  # Only apply s() if the column is numeric AND has enough unique values
  terms <- sapply(names(X_df), function(var) {
    if (is.numeric(X_df[[var]]) && length(unique(X_df[[var]])) > cts.num) {
      return(paste0("s(", var, ")"))
    } else {
      return(var)
    }
  })

  form_str <- paste("time ~", paste(terms, collapse = " + "))
  gam.formula <- as.formula(form_str)

  # 2. Prepare Data
  # mgcv::cox.ph requires 'time' in the data, and the event indicator passed as 'weights'
  dat <- data.frame(time = time, X_df)

  # 3. Fit GAM
  fit <- mgcv::gam(
    gam.formula,
    family = mgcv::cox.ph(),
    data = dat,
    weights = event, # Critical for mgcv cox models
    ...
  )

  # 4. Predict on newdata
  # mgcv natively returns S(t) if type="response" and newdata contains a 'time' column.
  n_new <- nrow(newdata_df)
  n_t   <- length(new.times)

  # Duplicate each row of newdata for every time point: P1, P1, P1, P2, P2, P2...
  long_newdata <- newdata_df[rep(1:n_new, each = n_t), , drop = FALSE]
  # Add the matching times: T1, T2, T3, T1, T2, T3...
  long_newdata$time <- rep(new.times, times = n_new)

  pred_vec <- predict(fit, newdata = long_newdata, type = "response")

  # Fill matrix by row to match the P1(T1,T2,T3) expansion structure
  pred <- matrix(pred_vec, nrow = n_new, ncol = n_t, byrow = TRUE)

  # 5. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  fit_obj <- list(object = fit, times = new.times)
  class(fit_obj) <- c("surv.gam")

  list(pred = pred, fit = fit_obj)
}



#' Prediction function for GAM wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.gam} object.
#'
#' @param object Fitted \code{surv.gam} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.gam(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     cts.num = 5
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.gam <- function(object, newdata, new.times, ...) {

  newdata_df <- as.data.frame(newdata)
  n_new <- nrow(newdata_df)
  n_t   <- length(new.times)

  # 1. Expand newdata to match mgcv's required format
  long_newdata <- newdata_df[rep(1:n_new, each = n_t), , drop = FALSE]
  long_newdata$time <- rep(new.times, times = n_new)

  # 2. Predict S(t)
  pred_vec <- predict(object$object, newdata = long_newdata, type = "response")
  pred <- matrix(pred_vec, nrow = n_new, ncol = n_t, byrow = TRUE)

  # 3. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}






#' Wrapper for standard Cox Proportional Hazards
#'
#' Final Production Wrapper for CoxPH.
#' Uses partial maximum likelihood and the Breslow estimator.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param ... Additional arguments passed to \code{\link[survival]{coxph}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.coxph(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' dim(fit$pred)
#' @export
surv.coxph <- function(time, event, X, newdata, new.times, obsWeights, id, ...) {
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Safety Clamp for Time
  time <- pmax(time, 1e-5)

  # 2. Clean Data Frame (Prevents '.' from swallowing external variables)
  dat <- data.frame(time = time, event = event, as.data.frame(X))

  # 3. Fit Cox Model safely passing 'id' only if it exists
  if(!missing(id) && !is.null(id)) {
    fit.coxph <- survival::coxph(
      survival::Surv(time, event) ~ .,
      data = dat,
      weights = obsWeights,
      id = id,
      x = TRUE,
      ...
    )
  } else {
    fit.coxph <- survival::coxph(
      survival::Surv(time, event) ~ .,
      data = dat,
      weights = obsWeights,
      x = TRUE,
      ...
    )
  }

  # 4. Get Baseline Hazard (Uncentered, evaluated at X=0)
  bh_clean <- survival::basehaz(fit.coxph, centered = FALSE)

  # 5. Fast Step Function Interpolation to new.times
  bh_fun <- stats::stepfun(bh_clean$time, c(0, bh_clean$hazard), right = FALSE)
  bh_interp <- bh_fun(new.times)

  # 6. Get Risk Scores on newdata (Must align with X=0 reference!)
  lp_new <- predict(
    fit.coxph,
    newdata = as.data.frame(newdata),
    type = "lp",
    reference = "zero"
  )

  # 7. Calculate Survival Probabilities
  pred <- outer(lp_new, bh_interp, function(eta, H0) exp(-exp(eta) * H0))

  # 8. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  # Save the raw basehaz table so predict.surv.coxph can evaluate any future time points
  fit <- list(
    object = fit.coxph,
    basehaz = bh_clean$hazard,
    times = bh_clean$time
  )
  class(fit) <- c("surv.coxph")

  list(pred = pred, fit = fit)
}




#' Prediction function for Cox regression wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.coxph} object.
#'
#' @param object Fitted \code{surv.coxph} object.
#' @param newdata New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:30, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:5, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.coxph(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' pred <- predict(fit$fit, newdata = newX, new.times = times)
#' dim(pred)
#' @export
predict.surv.coxph <- function(object, newdata, new.times, ...) {

  # 1. Get Linear Predictor (Risk Score) relative to X=0
  lp <- predict(
    object$object,
    newdata = as.data.frame(newdata),
    type = "lp",
    reference = "zero"
  )

  # 2. Reconstruct Baseline Hazard Step Function from saved raw data
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  H0 <- bh_fun(new.times)

  # 3. Calculate Survival S(t)
  pred <- outer(lp, H0, function(eta, h) exp(-exp(eta) * h))

  # 4. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}






#' Wrapper function for Component-Wise Boosting (CoxBoost)
#'
#' Final Production Wrapper for CoxBoost (Tunable & Robust).
#' Estimates a Cox model via component-wise likelihood based boosting.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights (Note: CoxBoost does not natively support weights, so these are ignored).
#' @param id Optional cluster/individual ID indicator.
#' @param stepno Number of boosting steps (default: 100).
#' @param penalty Penalty value for the update (default: 100).
#' @param ... Additional arguments passed to \code{\link[CoxBoost]{CoxBoost}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (requireNamespace("CoxBoost", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.coxboost(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     stepno = 10,
#'     penalty = 50
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.coxboost <- function(time, event, X, newdata, new.times, obsWeights, id,
                          stepno = 100, penalty = 100, ...) {

  requireNamespace("CoxBoost", quietly = TRUE)

  # 1. Prepare Matrices & Align Columns
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(colnames(X_mat), colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, colnames(X_mat), drop = FALSE]

  # 2. Fit CoxBoost
  # Note: CoxBoost automatically standardizes and selects variables step-by-step
  fit <- CoxBoost::CoxBoost(
    time = time,
    status = event,
    x = X_mat,
    stepno = stepno,
    penalty = penalty,
    ...
  )

  # 3. Predict Risk Scores (Linear Predictors)
  # CoxBoost returns a matrix for 'lp' if multiple steps are requested,
  # but by default returns the final step. We ensure it's a vector.
  lp_train <- as.numeric(predict(fit, newdata = X_mat, type = "lp"))

  # 4. Center Training Scores (To stabilize hazard calculation)
  tr_mean <- mean(lp_train)
  lp_train_centered <- lp_train - tr_mean

  # 5. Baseline Hazard (Using our robust helper)
  # Ensure safe_breslow_step is loaded in your environment!
  bh <- safe_breslow_step(
    time = time,
    event = event,
    risk_score = lp_train_centered,
    new.times = new.times
  )

  # 6. Predict on New Data
  lp_new <- as.numeric(predict(fit, newdata = newdata_mat, type = "lp"))
  lp_new_centered <- lp_new - tr_mean

  # SAFETY CLAMP: Prevent Inf * 0 = NaN
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 7. Convert to Survival Probabilities
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  # Save the column names so the predict method knows exactly what to expect
  fit_obj <- list(
    object = fit,
    basehaz = bh,
    times = new.times,
    stats = list(mean = tr_mean, features = colnames(X_mat))
  )
  class(fit_obj) <- c("surv.coxboost")

  list(pred = pred, fit = fit_obj)
}


#' Prediction function for CoxBoost wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.coxboost} object.
#'
#' @param object Fitted \code{surv.coxboost} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("CoxBoost", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.coxboost(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     stepno = 10,
#'     penalty = 50
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.coxboost <- function(object, newdata, new.times, ...) {

  # 1. Faster/Safer Interpolation using stepfun
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Format Matrix and Align Columns
  expected_cols <- object$stats$features
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(expected_cols, colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, expected_cols, drop = FALSE]

  # 3. Predict Linear Predictor
  lp_new <- as.numeric(predict(object$object, newdata = newdata_mat, type = "lp"))

  # 4. Center and Clamp
  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 5. Survival Formula
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newdata), ncol = length(new.times))

  # 6. Safety Clamping
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}








#' Wrapper for BART (Bayesian Additive Regression Trees)
#'
#' Final Production Wrapper for BART (Tunable & Robust).
#' Uses the \code{\link[BART]{mc.surv.bart}} function.
#' Automatically reshapes the flat output vector into a survival matrix
#' and interpolates the predictions to the requested new.times.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights (Note: BART does not natively support weights).
#' @param id Optional cluster/individual ID indicator.
#' @param ntree Number of trees (default: 50).
#' @param ndpost Number of posterior draws (default: 1000).
#' @param nskip Number of burn-in draws (default: 250).
#' @param ... Additional arguments passed to \code{\link[BART]{mc.surv.bart}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @examples
#' if (.Platform$OS.type != "windows" &&
#'   requireNamespace("BART", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:20, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.bart(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     ntree = 3,
#'     ndpost = 5,
#'     nskip = 5
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.bart <- function(time, event, X, newdata = NULL, new.times, obsWeights = NULL, id = NULL,
                      ntree = 10, ndpost = 30, nskip = 10, ...) {

  if (.Platform$OS.type == "windows") {
    stop(
      "surv.bart() is not supported on Windows because ",
      "BART::mc.surv.bart() requires multicore forking."
    )
  }

  requireNamespace("BART", quietly = TRUE)

  # 0) Defaults (important for SL internals)
  if (is.null(newdata)) newdata <- X

  # 1) Time safety (avoid exact 0)
  time <- pmax(time, 1e-5)

  # 2) Design matrices + column alignment
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(colnames(X_mat), colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, colnames(X_mat), drop = FALSE]

  # 3) Fit + predict (BART predicts x.test during fitting)
  fit <- BART::mc.surv.bart(
    x.train  = X_mat,
    times    = time,
    delta    = event,
    x.test   = newdata_mat,
    ntree    = ntree,
    ndpost   = ndpost,
    nskip    = nskip,
    mc.cores = 1,
    seed     = 99,
    ...
  )

  # 4) Extract BART-native time grid
  unique_times <- fit$times
  unique_times <- sort(unique_times)
  unique_times <- unique_times[unique_times > 0]

  n_new <- nrow(newdata_mat)
  n_t   <- length(unique_times)

  # 5) Reshape survival output robustly
  surv_vec <- fit$surv.test.mean
  if (length(surv_vec) != n_new * n_t) {
    stop(
      "BART returned an unexpected prediction length: ",
      length(surv_vec), " (expected ", n_new * n_t, ")."
    )
  }

  surv_raw <- matrix(surv_vec, nrow = n_new, ncol = n_t, byrow = TRUE)

  # 6) Interpolate to requested new.times (step-function style)
  #    yleft=1 ensures S(t)=1 before first event time.
  pred <- t(apply(surv_raw, 1, function(y) {
    stats::approx(
      x = unique_times, y = y, xout = new.times,
      method = "constant", f = 0, rule = 2, ties = mean,
      yleft = 1
    )$y
  }))

  # 7) Safety clamps + monotonicity
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) pred <- t(apply(pred, 1, cummin))

  fit_obj <- list(
    object = fit,
    times  = unique_times,
    stats  = list(features = colnames(X_mat), n_t = n_t)
  )
  class(fit_obj) <- c("surv.bart")

  list(pred = pred, fit = fit_obj)
}



#' Prediction function for BART
#'
#' Obtains predicted survivals from a fitted \code{surv.bart} object.
#' Manually expands the time grid to bypass the 10-column expectation.
#'
#' @param object Fitted \code{surv.bart} object.
#' @param newdata New covariate data.frame to predict on.
#' @param new.times Times to predict.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (.Platform$OS.type != "windows" &&
#'   requireNamespace("BART", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:20, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.bart(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     ntree = 3,
#'     ndpost = 5,
#'     nskip = 5
#'   )
#'
#'   pred <- fit$pred
#'   dim(pred)
#' }
#' @export
predict.surv.bart <- function(object, newdata, new.times, ...) {

  # 1) Align features
  expected_cols <- object$stats$features
  newdata_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newdata))

  missing_cols <- setdiff(expected_cols, colnames(newdata_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newdata_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newdata_mat <- cbind(newdata_mat, pad_mat)
  }
  newdata_mat <- newdata_mat[, expected_cols, drop = FALSE]

  # 2) Predict survival on the BART-native time grid
  #    (This is the least fragile approach; avoid manual "(t, X)" expansion.)
  p_obj <- tryCatch(
    {
      # Many BART objects support x.test style
      predict(object$object, x.test = newdata_mat)
    },
    error = function(e) NULL
  )

  if (is.null(p_obj)) {
    stop(
      "BART survival predictions on new data are not supported for this fitted object. ",
      "In SuperSurv, the BART wrapper is intended primarily for CV-time predictions ",
      "(via x.test during model fitting). If you need out-of-sample predictions, ",
      "consider using other learners or refitting BART with x.test=newdata inside the wrapper."
    )
  }

  # 3) Extract time grid + survival
  # Try common fields (these vary by version)
  if (!is.null(p_obj$times)) {
    bart_times <- p_obj$times
  } else if (!is.null(object$times)) {
    bart_times <- object$times
  } else {
    stop("Cannot find BART prediction time grid.")
  }

  # Survival output candidates
  surv_raw <- NULL
  if (!is.null(p_obj$surv.test.mean)) {
    # common
    surv_raw <- p_obj$surv.test.mean
  } else if (!is.null(p_obj$surv.mean)) {
    surv_raw <- p_obj$surv.mean
  } else if (!is.null(p_obj$surv)) {
    surv_raw <- p_obj$surv
  }

  if (is.null(surv_raw)) stop("Cannot find survival predictions in BART predict output.")

  # Ensure matrix shape: n_new x n_time
  n_new <- nrow(newdata_mat)
  n_t   <- length(bart_times)
  if (is.null(dim(surv_raw))) {
    surv_raw <- matrix(surv_raw, nrow = n_new, ncol = n_t, byrow = TRUE)
  }

  # 4) Interpolate to requested new.times (step function, survival monotone)
  pred <- t(apply(surv_raw, 1, function(y) {
    stats::stepfun(bart_times, c(1, y), right = FALSE)(new.times)
  }))

  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) pred <- t(apply(pred, 1, cummin))

  pred
}









#' Wrapper for AORSF (Oblique Random Survival Forest)
#'
#' Final Production Wrapper for AORSF (Tunable & Robust).
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newdata Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param n_tree Number of trees to grow (default: 500).
#' @param leaf_min_events Minimum number of events in a leaf node (default: 5).
#' @param mtry Number of predictors evaluated at each node.
#' @param ... Additional arguments passed to \code{\link[aorsf]{orsf}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }

#' @examples
#' if (requireNamespace("aorsf", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.aorsf(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     n_tree = 10,
#'     leaf_min_events = 2
#'   )
#'
#'   dim(fit$pred)
#' }
#' @export
surv.aorsf <- function(time, event, X, newdata, new.times, obsWeights, id,
                       n_tree = 500, leaf_min_events = 5, mtry = NULL, ...) {

  requireNamespace("aorsf", quietly = TRUE)

  # 1. Prepare Data
  dat <- data.frame(time = time, status = event, X)

  # Handle weights safely
  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 2. Fit Model
  # We construct a list of arguments for aorsf::orsf
  # We use this structure so that if mtry is NULL, aorsf uses its native default
  args_list <- list(
    data = dat,
    formula = time + status ~ .,
    weights = obsWeights,
    n_tree = n_tree,
    leaf_min_events = leaf_min_events,
    oobag_pred_type = "surv"  # Forced to 'surv' for SuperLearner compatibility
  )

  # Only add mtry if it was explicitly provided, otherwise let aorsf decide
  if (!is.null(mtry)) {
    args_list$mtry <- mtry
  }

  # Append any extra '...' arguments
  args_list <- c(args_list, list(...))

  # Execute the call safely
  fit <- do.call(aorsf::orsf, args_list)

  # 3. Predict Directly at new.times
  preds_surv <- predict(
    fit,
    new_data = newdata,
    pred_horizon = new.times,
    pred_type = "surv"
  )

  # 4. Format Output & Safety (Monotonicity)
  pred_mat <- as.matrix(preds_surv)
  if (ncol(pred_mat) > 1) {
    pred_mat <- t(apply(pred_mat, 1, cummin))
  }

  fit_obj <- list(object = fit)
  class(fit_obj) <- c("surv.aorsf")

  list(pred = pred_mat, fit = fit_obj)
}


#' Prediction function for AORSF
#'
#' Obtains predicted survivals from a fitted \code{surv.aorsf} object.
#' Uses the native \code{aorsf} prediction engine to calculate survival directly at requested times.
#'
#' @param object Fitted \code{surv.aorsf} object.
#' @param newdata New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newdata} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @examples
#' if (requireNamespace("aorsf", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:30, ]
#'   x_cols <- grep("^x", names(dat))[1:3]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:5, , drop = FALSE]
#'   times <- seq(50, 150, by = 50)
#'
#'   fit <- surv.aorsf(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = newX,
#'     new.times = times,
#'     obsWeights = rep(1, nrow(dat)),
#'     id = NULL,
#'     n_tree = 10,
#'     leaf_min_events = 2
#'   )
#'
#'   pred <- predict(fit$fit, newdata = newX, new.times = times)
#'   dim(pred)
#' }
#' @export
predict.surv.aorsf <- function(object, newdata, new.times, ...) {

  # 1. Predict directly at the requested new.times
  # object$object is the actual 'orsf' model
  preds_surv <- predict(
    object$object,
    new_data = newdata,
    pred_horizon = new.times,
    pred_type = "surv"
  )

  # 2. Format as Matrix
  out <- as.matrix(preds_surv)

  # 3. Safety Check: Dimensions
  # If new.times was a single number, ensure it's a column matrix
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = nrow(newdata), ncol = length(new.times))
  }

  out[out < 0] <- 0
  out[out > 1] <- 1

  # 4. Safety Check: Monotonicity
  # Ensure probability never increases over time (S(t) must go down)
  if (ncol(out) > 1) {
    out <- t(apply(out, 1, cummin))
  }

  return(out)
}










