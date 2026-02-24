
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
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.xgboost <- function(time, event, X, newX, new.times, obsWeights, id,
                         nrounds = 1000, early_stopping_rounds = 10,
                         eta = 0.05, max_depth = 2, min_child_weight = 5,
                         lambda = 10, subsample = 0.7, ...) {

  requireNamespace("xgboost", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  # 1. Prepare Matrices
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  # CRITICAL FIX: Matrix Alignment for CV
  # If newX is missing a factor level present in X, we must pad it with zeros
  missing_cols <- setdiff(colnames(X_mat), colnames(newX_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newX_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newX_mat <- cbind(newX_mat, pad_mat)
  }
  # Ensure columns are in the exact same order
  newX_mat <- newX_mat[, colnames(X_mat), drop = FALSE]

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
  lp_new <- predict(fit, newdata = newX_mat)
  lp_new_centered <- lp_new - tr_mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  # 7. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(object = fit, basehaz = bh, times = new.times, stats = list(mean = tr_mean))
  class(fit_obj) <- c("surv.xgboost")

  list(pred = pred, fit = fit_obj)
}





#' Prediction function for XGBoost wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.xgboost} object.
#'
#' @param object Fitted \code{surv.xgboost} object.
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.xgboost <- function(object, newX, new.times, ...) {

  # 1. Faster/Safer Interpolation using stepfun
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Format Matrix and Align Columns
  # Must use the exact features the XGBoost model expects
  expected_cols <- object$object$feature_names
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  missing_cols <- setdiff(expected_cols, colnames(newX_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newX_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newX_mat <- cbind(newX_mat, pad_mat)
  }
  newX_mat <- newX_mat[, expected_cols, drop = FALSE]

  # 3. Predict Linear Predictor
  lp_new <- predict(object$object, newdata = newX_mat)

  # 4. Center and Clamp
  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 5. Survival Formula
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  # 6. Safety Clamping
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}




