
#' Universal Parametric Survival Wrapper
#'
#' Final Production Wrapper for AFT Models (Weibull, Exponential, LogNormal, LogLogistic).
#' Replaces individual wrappers with one robust, vectorized function.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param dist Distribution for the AFT model (default: "weibull").
#' @param ... Additional arguments passed to \code{\link[survival]{survreg}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @export
surv.parametric <- function(time, event, X, newX, new.times, obsWeights, id, dist = "weibull", ...) {
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Handle Time <= 0 issue for AFT models (log(0) = -Inf)
  time <- pmax(time, 1e-5)

  # 2. Fit Model
  # CRITICAL FIX: Use data.frame() to prevent cbind from destroying factor variables
  dat <- data.frame(time = time, event = event, as.data.frame(X))

  fit <- survival::survreg(
    survival::Surv(time, event) ~ .,
    data = dat,
    weights = obsWeights,
    dist = dist,
    ...
  )

  # 3. Fast Prediction (Closed Form)
  # predict.survreg uses type = "linear" to return the linear predictor
  lp <- predict(fit, newdata = as.data.frame(newX), type = "linear")
  scale <- fit$scale

  # 4. Vectorized Survival Calculation
  if (dist %in% c("exponential", "weibull")) {
    pred <- outer(lp, new.times, function(eta, t) {
      exp( - (t / exp(eta))^(1/scale) )
    })

  } else if (dist == "lognormal") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 - stats::pnorm((log(t) - eta) / scale)
    })

  } else if (dist == "loglogistic") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 / (1 + (t / exp(eta))^(1/scale))
    })

  } else {
    stop("Distribution not supported in fast wrapper")
  }

  # 5. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1

  fit_obj <- list(object = fit, dist = dist)
  class(fit_obj) <- c("surv.parametric")

  list(pred = pred, fit = fit_obj)
}

# --- Specific Wrappers Pointing to the Universal One ---

#' Parametric Survival Prediction Wrapper (Exponential)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @export
surv.exponential <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newX = newX,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "exponential", ...)
}

#' Parametric Survival Prediction Wrapper (Log-Logistic)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @export
surv.loglogistic <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newX = newX,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "loglogistic", ...)
}

#' Parametric Survival Prediction Wrapper (Log-Normal)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @export
surv.lognormal <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newX = newX,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "lognormal", ...)
}

#' Parametric Survival Prediction Wrapper (Weibull)
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Cluster identification variable.
#' @param ... Additional ignored arguments.
#'
#' @return A list containing the fitted model and predictions.
#' @export
surv.weibull <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  surv.parametric(time = time, event = event, X = X, newX = newX,
                  new.times = new.times, obsWeights = obsWeights, id = id,
                  dist = "weibull", ...)
}






#' Prediction function for Universal Parametric Wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.parametric} object
#' using exact closed-form equations.
#'
#' @param object Fitted \code{surv.parametric} object.
#' @param newX New covariate data.frame.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.parametric <- function(object, newX, new.times, ...) {

  fit <- object$object
  dist <- object$dist
  scale <- fit$scale

  # Extract linear predictor
  lp <- predict(fit, newdata = as.data.frame(newX), type = "linear")

  # Vectorized Survival Calculation
  if (dist %in% c("exponential", "weibull")) {
    pred <- outer(lp, new.times, function(eta, t) {
      exp( - (t / exp(eta))^(1/scale) )
    })
  } else if (dist == "lognormal") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 - stats::pnorm((log(t) - eta) / scale)
    })
  } else if (dist == "loglogistic") {
    pred <- outer(lp, new.times, function(eta, t) {
      1 / (1 + (t / exp(eta))^(1/scale))
    })
  }

  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  return(pred)
}

