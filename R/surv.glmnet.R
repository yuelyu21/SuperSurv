

############################################################
## GLMNET (Lasso/ElasticNet) — wrapper + predict method
## Pattern: Risk Score + Cox Offset Calibration
############################################################


#' Wrapper function for Penalized Cox Regression (GLMNET)
#'
#' Final Production Wrapper for GLMNET (Tunable & Robust).
#' Estimates a penalized Cox model (Lasso, Ridge, or Elastic Net) with automatic lambda selection.
#' Uses the Breslow estimator with a step-function approach for the baseline hazard.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.glmnet <- function(time, event, X, newX, new.times, obsWeights, id,
                        alpha = 1, nfolds = 10, ...) {

  requireNamespace("glmnet", quietly = TRUE)

  # We clamp the minimum time to a tiny positive number.
  time <- pmax(time, 1e-5)

  # 1. Convert factors to dummy variables (glmnet requires matrices)
  X_mat <- stats::model.matrix(~ . - 1, data = X)
  newX_mat <- stats::model.matrix(~ . - 1, data = newX)

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

  # 3. Calibrate Baseline Hazard (Breslow Estimator)
  lp_train <- as.numeric(predict(fit, newx = X_mat, s = "lambda.min", type = "link"))
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

  # 4. Predict on newX
  lp_new <- as.numeric(predict(fit, newx = newX_mat, s = "lambda.min", type = "link"))
  pred <- outer(lp_new, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  # 5. Safety Clamp
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(object = fit, basehaz = bh, times = new.times)
  class(fit_obj) <- c("surv.glmnet")

  list(pred = pred, fit = fit_obj)
}




#' Prediction function for GLMNET wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.glmnet} object.
#'
#' @param object Fitted \code{surv.glmnet} object.
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.glmnet <- function(object, newX, new.times, ...) {

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

  # 2. Convert newX to matrix safely
  # model.matrix removes NAs, so we ensure the dataframe is intact
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  # 3. Predict Linear Predictor
  lp_new <- as.numeric(predict(
    object$object,
    newx = newX_mat,
    s = "lambda.min",
    type = "link"
  ))

  # 4. Convert LP to Survival Probability Matrix
  pred <- outer(lp_new, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  # 5. Safety Clamp
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}




