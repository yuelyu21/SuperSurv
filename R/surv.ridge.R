#' Wrapper for Ridge Regression (Penalized Cox)
#'
#' Final Production Wrapper for Ridge Regression (Tunable & Robust).
#' Estimates a penalized Cox model using a pure Ridge penalty (alpha = 0).
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.ridge <- function(time, event, X, newX, new.times, obsWeights = NULL, id = NULL, nfolds = 10, ...) {

  requireNamespace("glmnet", quietly = TRUE)

  # We clamp the minimum time to a tiny positive number.
  time <- pmax(time, 1e-5)

  if (is.null(obsWeights)) obsWeights <- rep(1, length(time))

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

  # 5. Predict on newX
  lp_new <- as.numeric(predict(fit, newx = newX_mat, s = "lambda.min", type = "link"))
  lp_new_centered <- lp_new - tr_mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.ridge <- function(object, newX, new.times, ...) {

  # 1. Faster/Safer Interpolation
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
  lp_new <- as.numeric(predict(object$object, newx = newX_mat, s = "lambda.min", type = "link"))

  # 4. Center and Clamp
  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 5. Survival Formula
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  pred[pred < 0] <- 0
  pred[pred > 1] <- 1

  return(pred)
}

