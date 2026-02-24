#' Wrapper for standard Cox Proportional Hazards
#'
#' Final Production Wrapper for CoxPH.
#' Uses partial maximum likelihood and the Breslow estimator.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Times at which to obtain predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Optional cluster/individual ID indicator.
#' @param ... Additional arguments passed to \code{\link[survival]{coxph}}.
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: The fitted model object (e.g., the raw \code{coxph} or \code{xgb.Booster} object).
#'     If the model fails to fit, this may be an object of class \code{try-error}.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions
#'     evaluated at the specified \code{new.times} grid.
#' }
#' @export
surv.coxph <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  requireNamespace("survival", quietly = TRUE)

  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))

  # 1. Safety Clamp for Time
  time <- pmax(time, 1e-5)

  # 2. Clean Data Frame (Prevents '.' from swallowing external variables)
  dat <- data.frame(time = time, event = event, as.data.frame(X))

  # 3. Fit Cox Model safely passing 'id' only if it exists
  if(!missing(id) && !is.null(id)) {
    fit.coxph <- survival::coxph(
      survival::Surv(time, event) ~ .,
      data = dat,
      weights = obsWeights,
      id = id,
      x = TRUE,
      ...
    )
  } else {
    fit.coxph <- survival::coxph(
      survival::Surv(time, event) ~ .,
      data = dat,
      weights = obsWeights,
      x = TRUE,
      ...
    )
  }

  # 4. Get Baseline Hazard (Uncentered, evaluated at X=0)
  bh_clean <- survival::basehaz(fit.coxph, centered = FALSE)

  # 5. Fast Step Function Interpolation to new.times
  bh_fun <- stats::stepfun(bh_clean$time, c(0, bh_clean$hazard), right = FALSE)
  bh_interp <- bh_fun(new.times)

  # 6. Get Risk Scores on newX (Must align with X=0 reference!)
  lp_new <- predict(
    fit.coxph,
    newdata = as.data.frame(newX),
    type = "lp",
    reference = "zero"
  )

  # 7. Calculate Survival Probabilities
  pred <- outer(lp_new, bh_interp, function(eta, H0) exp(-exp(eta) * H0))

  # 8. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  # Save the raw basehaz table so predict.surv.coxph can evaluate any future time points
  fit <- list(
    object = fit.coxph,
    basehaz = bh_clean$hazard,
    times = bh_clean$time
  )
  class(fit) <- c("surv.coxph")

  list(pred = pred, fit = fit)
}




#' Prediction function for Cox regression wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.coxph} object.
#'
#' @param object Fitted \code{surv.coxph} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain predicted survivals.
#' @param ... Additional ignored arguments.
#' @return A numeric matrix of predicted survival probabilities, where rows correspond
#'   to the observations in \code{newX} and columns correspond to the evaluation
#'   times in \code{new.times}.
#' @export
predict.surv.coxph <- function(object, newX, new.times, ...) {

  # 1. Get Linear Predictor (Risk Score) relative to X=0
  lp <- predict(
    object$object,
    newdata = as.data.frame(newX),
    type = "lp",
    reference = "zero"
  )

  # 2. Reconstruct Baseline Hazard Step Function from saved raw data
  bh_fun <- stats::stepfun(object$times, c(0, object$basehaz), right = FALSE)
  H0 <- bh_fun(new.times)

  # 3. Calculate Survival S(t)
  pred <- outer(lp, H0, function(eta, h) exp(-exp(eta) * h))

  # 4. Safety Clamps
  pred[pred < 0] <- 0; pred[pred > 1] <- 1
  if (ncol(pred) > 1) {
    pred <- t(apply(pred, 1, cummin))
  }

  return(pred)
}




