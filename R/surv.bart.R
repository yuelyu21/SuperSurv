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
#' @param newX Test covariate data.frame to use for prediction.
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
#' @export
surv.bart <- function(time, event, X, newX, new.times, obsWeights, id,
                      ntree = 30, ndpost = 200, nskip = 100, ...) {

  requireNamespace("BART", quietly = TRUE)

  # ==========================================
  # THE BULLETPROOF TIME BLOCK
  # 1. Round first (Optional: speeds up BART heavily)
  time <- round(time, 1)

  # 2. Clamp SECOND (Ensures no exact zeroes ever exist)
  time <- pmax(time, 1e-5)
  # 2. Prepare Matrices & Align Columns
  X_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(X))
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  missing_cols <- setdiff(colnames(X_mat), colnames(newX_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newX_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newX_mat <- cbind(newX_mat, pad_mat)
  }
  newX_mat <- newX_mat[, colnames(X_mat), drop = FALSE]

  # 3. Fit & Predict (BART runs newX during fitting)
  # mc.cores = 1 is standard to prevent nested parallelization crashes in SuperLearner
  fit <- BART::mc.surv.bart(
    x.train = X_mat,
    times   = time,
    delta   = event,
    x.test  = newX_mat,
    ntree   = ntree,
    ndpost  = ndpost,
    nskip   = nskip,
    mc.cores = 1,
    seed    = 99,
    ...
  )

  # 4. Reshape Flat Output
  unique_times <- fit$times
  n_new <- nrow(newX_mat)
  n_t   <- length(unique_times)

  surv_raw <- matrix(
    fit$surv.test.mean,
    nrow = n_new,
    ncol = n_t,
    byrow = TRUE
  )

  # 5. Vectorized Interpolation to new.times
  # stepfun ensures survival at t=0 is 1.0
  pred <- t(apply(surv_raw, 1, function(y) {
    stats::stepfun(unique_times, c(1, y), right = FALSE)(new.times)
  }))

  # 6. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  fit_obj <- list(
    object = fit,
    times = unique_times,
    stats = list(features = colnames(X_mat))
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
#' @param newX New covariate data.frame to predict on.
#' @param new.times Times to predict.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.bart <- function(object, newX, new.times, ...) {

  # 1. Format Matrix and Align Columns (9 columns for METABRIC)
  expected_cols <- object$stats$features
  newX_mat <- stats::model.matrix(~ . - 1, data = as.data.frame(newX))

  missing_cols <- setdiff(expected_cols, colnames(newX_mat))
  if (length(missing_cols) > 0) {
    pad_mat <- matrix(0, nrow = nrow(newX_mat), ncol = length(missing_cols))
    colnames(pad_mat) <- missing_cols
    newX_mat <- cbind(newX_mat, pad_mat)
  }
  newX_mat <- newX_mat[, expected_cols, drop = FALSE]

  # 2. Extract unique times from the trained object
  unique_times <- object$times
  n_new <- nrow(newX_mat)
  n_t   <- length(unique_times)

  # 3. Create Expanded Matrix
  # Repeat each patient row for every unique time
  expanded_newX <- newX_mat[rep(1:n_new, each = n_t), , drop = FALSE]

  # Add the time column 't' at the very beginning (Making it 10 columns!)
  t_col <- rep(unique_times, times = n_new)
  expanded_newX <- cbind(t = t_col, expanded_newX)

  # 4. Predict discrete hazards using the main object directly
  pbart_pred <- predict(object$object, newdata = expanded_newX)

  # 5. Extract the mean hazard across all posterior draws
  # Depending on the BART version, it returns prob.test or prob.test.mean
  if (!is.null(pbart_pred$prob.test)) {
    hazards <- colMeans(pbart_pred$prob.test)
  } else {
    hazards <- pbart_pred$prob.test.mean
  }

  # 6. Reshape and Convert Hazards to Survival S(t)
  haz_mat <- matrix(hazards, nrow = n_new, ncol = n_t, byrow = TRUE)

  # S(t) = cumulative product of (1 - discrete hazard)
  surv_raw <- t(apply(haz_mat, 1, function(h) cumprod(1 - h)))

  # 7. Vectorized Interpolation to new.times
  pred <- t(apply(surv_raw, 1, function(y) {
    stats::stepfun(unique_times, c(1, y), right = FALSE)(new.times)
  }))

  # 8. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}


