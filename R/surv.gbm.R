############################################################
## GBM (coxph) — wrapper + predict method (robust)
############################################################

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
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.gbm <- function(time, event, X, newX, new.times, obsWeights, id,
                     n.trees = 1000, interaction.depth = 2, shrinkage = 0.01,
                     cv.folds = 5, n.minobsinnode = 10, ...) {

  requireNamespace("gbm", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  # 1. CRITICAL FIX: Prevent C++ Crashes (Data Types)
  if(is.matrix(X)) X <- as.data.frame(X)
  if(is.matrix(newX)) newX <- as.data.frame(newX)

  # Force characters to factors
  X[] <- lapply(X, function(x) if(is.character(x)) as.factor(x) else x)
  newX[] <- lapply(newX, function(x) if(is.character(x)) as.factor(x) else x)

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
  lp_new <- predict(fit, newdata = newX, n.trees = best.iter, type = "link")

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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.gbm <- function(object, newX, new.times, ...) {

  # 1. Safer Interpolation using stepfun
  # Prepend 0 to handle times before the first event safely
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Linear Predictor
  lp_new <- predict(
    object$object,
    newdata = as.data.frame(newX),
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








