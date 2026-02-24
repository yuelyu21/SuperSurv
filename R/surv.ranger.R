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
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.ranger <- function(time, event, X, newX, new.times, obsWeights, id,
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

  # 4. Predict on newX (Training Grid)
  p_obj <- predict(fit, data = as.data.frame(newX))
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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.ranger <- function(object, newX, new.times, ...) {

  # 1. Predict using saved model
  p_obj <- predict(object$object, data = as.data.frame(newX))
  surv_probs  <- p_obj$survival
  train_times <- p_obj$unique.death.times

  # 2. Vectorized Interpolation
  pred <- t(apply(surv_probs, 1, function(y) {
    stats::approx(x = train_times, y = y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  # 3. Dimensions Safety Check
  if (is.null(dim(pred))) {
    pred <- matrix(pred, nrow = nrow(newX), ncol = length(new.times))
  }

  # 4. Safety Clamp and Monotonicity
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}



