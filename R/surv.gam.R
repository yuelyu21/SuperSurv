#' Wrapper for Generalized Additive Cox Regression (GAM)
#'
#' Final Production Wrapper for GAM (Tunable & Robust).
#' Uses \code{\link[mgcv]{gam}} to fit an additive combination of smooth
#' and linear functions.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.gam <- function(time, event, X, newX, new.times, obsWeights, id, cts.num = 5, ...) {

  requireNamespace("mgcv", quietly = TRUE)

  # 1. Formula Construction (Safe Spline Assignment)
  X_df <- as.data.frame(X)
  newX_df <- as.data.frame(newX)

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

  # 4. Predict on newX
  # mgcv natively returns S(t) if type="response" and newdata contains a 'time' column.
  n_new <- nrow(newX_df)
  n_t   <- length(new.times)

  # Duplicate each row of newX for every time point: P1, P1, P1, P2, P2, P2...
  long_newX <- newX_df[rep(1:n_new, each = n_t), , drop = FALSE]
  # Add the matching times: T1, T2, T3, T1, T2, T3...
  long_newX$time <- rep(new.times, times = n_new)

  pred_vec <- predict(fit, newdata = long_newX, type = "response")

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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.gam <- function(object, newX, new.times, ...) {

  newX_df <- as.data.frame(newX)
  n_new <- nrow(newX_df)
  n_t   <- length(new.times)

  # 1. Expand newX to match mgcv's required format
  long_newX <- newX_df[rep(1:n_new, each = n_t), , drop = FALSE]
  long_newX$time <- rep(new.times, times = n_new)

  # 2. Predict S(t)
  pred_vec <- predict(object$object, newdata = long_newX, type = "response")
  pred <- matrix(pred_vec, nrow = n_new, ncol = n_t, byrow = TRUE)

  # 3. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}



