#' Wrapper function for Component-Wise Boosting (CoxBoost)
#'
#' Final Production Wrapper for CoxBoost (Tunable & Robust).
#' Estimates a Cox model via component-wise likelihood based boosting.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.coxboost <- function(time, event, X, newX, new.times, obsWeights, id,
                          stepno = 100, penalty = 100, ...) {

  requireNamespace("CoxBoost", quietly = TRUE)

  # 1. Prepare Matrices & Align Columns
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  missing_cols <- setdiff(colnames(X_mat), colnames(newX_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newX_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newX_mat <- cbind(newX_mat, pad_mat)
  }
  newX_mat <- newX_mat[, colnames(X_mat), drop = FALSE]

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
  lp_new <- as.numeric(predict(fit, newdata = newX_mat, type = "lp"))
  lp_new_centered <- lp_new - tr_mean

  # SAFETY CLAMP: Prevent Inf * 0 = NaN
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 7. Convert to Survival Probabilities
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.coxboost <- function(object, newX, new.times, ...) {

  # 1. Faster/Safer Interpolation using stepfun
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Format Matrix and Align Columns
  expected_cols <- object$stats$features
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  missing_cols <- setdiff(expected_cols, colnames(newX_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newX_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newX_mat <- cbind(newX_mat, pad_mat)
  }
  newX_mat <- newX_mat[, expected_cols, drop = FALSE]

  # 3. Predict Linear Predictor
  lp_new <- as.numeric(predict(object$object, newdata = newX_mat, type = "lp"))

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




