#' Wrapper function for Random Survival Forests (RFSRC)
#'
#' Final Production Wrapper for RFSRC (Tunable & Robust).
#' Estimates a survival random forest using \code{\link[randomForestSRC]{rfsrc}}.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Currently ignored.
#' @param ntree Number of trees to grow (default: 1000).
#' @param nodesize Minimum number of deaths in terminal nodes (default: 15).
#' @param mtry Number of variables randomly selected as candidates for splitting a node.
#' @param ... Additional arguments passed to \code{\link[randomForestSRC]{rfsrc}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @export
surv.rfsrc <- function(time, event, X, newX, new.times, obsWeights = NULL, id = NULL,
                       ntree = 1000, nodesize = 15, mtry = NULL, ...) {

  requireNamespace("randomForestSRC", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  # 1. Handle Defaults
  if(is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 2. Prepare Data (RFSRC strictly needs data.frames)
  X_df <- as.data.frame(X)
  newX_df <- as.data.frame(newX)
  data <- data.frame(time = time, event = event, X_df)

  # 3. Setup Arguments
  args_list <- list(
    formula = survival::Surv(time, event) ~ .,
    data = data,
    case.wt = obsWeights,
    ntree = ntree,
    nodesize = nodesize,
    importance = "none", # Faster fitting
    ...
  )

  # Only add mtry if explicitly provided, otherwise let rfsrc calculate default
  if (!is.null(mtry)) args_list$mtry <- mtry

  # 4. Fit Model
  fit.rfsrc <- do.call(randomForestSRC::rfsrc, args_list)

  # 5. Predict on newX (Training Grid)
  survs <- predict(fit.rfsrc, newdata = newX_df, importance = "none")$survival
  train_times <- fit.rfsrc$time.interest

  # 6. Interpolate to new.times
  pred <- t(apply(survs, 1, function(y) {
    stats::approx(train_times, y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  # 7. Safety: Monotonicity and Clamping
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  fit <- list(object = fit.rfsrc)
  class(fit) <- c("surv.rfsrc")
  list(pred = pred, fit = fit)
}



#' Prediction function for survival random forest (RFSRC)
#'
#' Obtains predicted survivals from a fitted \code{surv.rfsrc} object.
#'
#' @param object Fitted \code{surv.rfsrc} object.
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.rfsrc <- function(object, newX, new.times, ...) {

  fit <- object$object

  # 1. Predict on native time grid
  survs <- predict(fit, newdata = as.data.frame(newX), importance = "none")$survival
  train_times <- fit$time.interest

  # 2. Vectorized interpolation to new.times
  pred <- t(apply(survs, 1, function(y) {
    stats::approx(train_times, y, xout = new.times, method = "constant", rule = 2, ties = mean)$y
  }))

  # 3. Safety Check: Dimensions
  if (is.null(dim(pred))) {
    pred <- matrix(pred, nrow = nrow(newX), ncol = length(new.times))
  }

  # 4. Safety Check: Clamp and Monotonicity
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}




