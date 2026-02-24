#' Kaplan-Meier Prediction Algorithm
#'
#' This prediction algorithm ignores all covariates and computes the marginal 
#' Kaplan-Meier survival estimator using the \code{\link[survival]{survfit}} function.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param X Training covariate data.frame (Ignored by KM).
#' @param newX Test covariate data.frame to use for prediction.
#' @param new.times Numeric vector of times at which to predict survival.
#' @param obsWeights Numeric vector of observation weights.
#' @param id Optional vector indicating subject/cluster identities.
#' @param ... Additional ignored arguments.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: A list containing the fitted \code{\link[survival]{survfit}} object.
#'   \item \code{pred}: A numeric matrix of cross-validated survival predictions 
#'     evaluated at \code{new.times}.
#' }
#' @export
surv.km <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  
  if(missing(obsWeights) || is.null(obsWeights)) obsWeights <- rep(1, length(time))
  
  # Fit the weighted KM curve
  fit.km <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    weights = obsWeights
  )
  
  # Efficient prediction using stepfun
  sfun <- stats::stepfun(fit.km$time, c(1, fit.km$surv), right = FALSE)
  surv_probs <- sfun(new.times)
  
  # Repeat for every patient (KM gives the exact same curve for everyone)
  pred <- matrix(surv_probs, 
                 nrow = nrow(newX), 
                 ncol = length(new.times), 
                 byrow = TRUE)
  
  fit <- list(object = fit.km)
  class(fit) <- c("surv.km")
  
  return(list(pred = pred, fit = fit))
}


#' Predict Method for Kaplan-Meier Wrapper
#'
#' Obtains predicted survivals from a fitted \code{surv.km} object.
#'
#' @param object A fitted object of class \code{surv.km}.
#' @param newX New covariate data.frame for which to obtain predictions (Ignored).
#' @param new.times Numeric vector of times at which to predict survival.
#' @param ... Additional ignored arguments.
#' 
#' @return A numeric matrix of predicted survival probabilities, where rows correspond 
#'   to the observations in \code{newX} and columns correspond to the evaluation times.
#' @export
#' @keywords internal
predict.surv.km <- function(object, newX, new.times, ...) {
  
  # Extract the original fitted KM object
  fit.km <- object$object
  
  # Reconstruct the step function
  sfun <- stats::stepfun(fit.km$time, c(1, fit.km$surv), right = FALSE)
  surv_probs <- sfun(new.times)
  
  # Replicate across all new patients
  pred_matrix <- matrix(surv_probs, 
                        nrow = nrow(newX), 
                        ncol = length(new.times), 
                        byrow = TRUE)
  
  return(pred_matrix)
}