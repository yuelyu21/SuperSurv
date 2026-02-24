#' Wrapper for Survival Regression Trees (rpart)
#'
#' Final Production Wrapper for single decision trees.
#' Uses \code{\link[rpart]{rpart}} with method="exp" and calculates
#' survival probabilities using the Breslow estimator.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param cp Complexity parameter (default: 0.01).
#' @param minsplit Minimum number of observations to attempt a split (default: 20).
#' @param maxdepth Maximum depth of any node of the final tree (default: 30).
#' @param ... Additional arguments passed to \code{\link[rpart]{rpart.control}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @export
surv.rpart <- function(time, event, X, newX, new.times, obsWeights, id,
                       cp = 0.01, minsplit = 20, maxdepth = 30, ...) {

  requireNamespace("rpart", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Safety Clamp for Time
  time <- pmax(time, 1e-5)

  # 2. Prepare Data
  X_df <- as.data.frame(X)
  newX_df <- as.data.frame(newX)
  dat <- data.frame(time = time, event = event, X_df)

  # 3. Fit the Survival Tree
  fit <- rpart::rpart(
    survival::Surv(time, event) ~ .,
    data = dat,
    weights = obsWeights,
    method = "exp", # Survival method for rpart
    control = rpart::rpart.control(
      cp = cp,
      minsplit = minsplit,
      maxdepth = maxdepth,
      ...
    )
  )

  # 4. Predict Expected Event Rates and Convert to Linear Predictor
  # If a leaf has 0 events, the rate is 0. log(0) = -Inf, which crashes Breslow.
  # We use pmax to safely clamp the minimum rate before taking the log.
  rate_train <- predict(fit, newdata = X_df)
  lp_train <- log(pmax(rate_train, 1e-10))

  # 5. Center Training Scores
  tr_mean <- mean(lp_train)
  lp_train_centered <- lp_train - tr_mean

  # 6. Baseline Hazard Calculation
  bh <- safe_breslow_step(
    time = time,
    event = event,
    risk_score = lp_train_centered,
    new.times = new.times
  )

  # 7. Predict on newX
  rate_new <- predict(fit, newdata = newX_df)
  lp_new <- log(pmax(rate_new, 1e-10))

  lp_new_centered <- lp_new - tr_mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 8. Calculate Survival S(t)
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX_df), ncol = length(new.times))

  # 9. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  fit_obj <- list(
    object = fit,
    basehaz = bh,
    times = new.times,
    stats = list(mean = tr_mean)
  )
  class(fit_obj) <- c("surv.rpart")

  list(pred = pred, fit = fit_obj)
}


#' Prediction function for rpart wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.rpart} object.
#'
#' @param object Fitted \code{surv.rpart} object.
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.rpart <- function(object, newX, new.times, ...) {

  # 1. Reconstruct Baseline Hazard
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  bh <- bh_fun(new.times)

  # 2. Predict Event Rate and Convert to Centered LP
  rate_new <- predict(object$object, newdata = as.data.frame(newX))
  lp_new <- log(pmax(rate_new, 1e-10))

  lp_new_centered <- lp_new - object$stats$mean
  lp_new_centered <- pmin(lp_new_centered, 700)

  # 3. Calculate Survival S(t)
  pred <- outer(lp_new_centered, bh, function(lp, h) exp(-exp(lp) * h))
  pred <- matrix(as.vector(pred), nrow = nrow(newX), ncol = length(new.times))

  # 4. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}

