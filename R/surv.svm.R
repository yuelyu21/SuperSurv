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
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.svm <- function(time, event, X, newX, new.times, obsWeights, id,
                     gamma.mu = 0.1, type = "vanbelle2",
                     kernel = "lin_kernel", opt.meth = "quadprog", ...) {

  requireNamespace("survivalsvm", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Prepare Data
  Y <- survival::Surv(time, event)
  df_train <- data.frame(Y = Y, X)

  # 2. Fit Survival SVM
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
  bh_df <- survival::basehaz(calib_fit, centered = TRUE)
  bh <- stats::approx(
    bh_df$time, bh_df$hazard, xout = new.times,
    method = "constant", f = 0, rule = 2, ties = mean
  )$y
  bh <- cummax(replace(bh, is.na(bh), 0))

  # 6. Predict and Calibrate on New Data
  p_te <- predict(fit, newdata = newX)
  lp_new_raw <- -1 * as.numeric(p_te$predicted)

  # Use R's native predict.coxph to get perfectly centered linear predictors
  lp_new_calibrated <- predict(
    calib_fit,
    newdata = data.frame(score = lp_new_raw),
    type = "lp"
  )

  # 7. Convert to Survival Probabilities
  pred <- outer(lp_new_calibrated, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.svm <- function(object, newX, new.times, ...) {

  # 1. Faster/Safer Interpolation using stepfun
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Get Raw Predictions from SVM
  p_te <- predict(object$object, newdata = newX)
  lp_new_raw <- -1 * as.numeric(p_te$predicted)

  # 3. Apply Calibration using the saved Cox model
  lp_new_calibrated <- predict(
    object$calibrator,
    newdata = data.frame(score = lp_new_raw),
    type = "lp"
  )

  # 4. Convert to Survival Probabilities S(t)
  pred <- outer(lp_new_calibrated, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  # 5. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  return(pred)
}









