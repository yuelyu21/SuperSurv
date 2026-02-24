############################################################
## AORSF — survSuperLearner wrapper + predict method
############################################################


#' Wrapper for AORSF (Oblique Random Survival Forest)
#'
#' Final Production Wrapper for AORSF (Tunable & Robust).
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
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
#' @export
surv.aorsf <- function(time, event, X, newX, new.times, obsWeights, id,
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
    new_data = newX,
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
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
#' @keywords internal
predict.surv.aorsf <- function(object, newX, new.times, ...) {

  # 1. Predict directly at the requested new.times
  # object$object is the actual 'orsf' model
  preds_surv <- predict(
    object$object,
    new_data = newX,
    pred_horizon = new.times,
    pred_type = "surv"
  )

  # 2. Format as Matrix
  out <- as.matrix(preds_surv)

  # 3. Safety Check: Dimensions
  # If new.times was a single number, ensure it's a column matrix
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = nrow(newX), ncol = length(new.times))
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



