#' Predict method for SuperSurv fits
#'
#' Obtains predicted survival probabilities from a fitted SuperSurv ensemble.
#'
#' @param object A fitted object of class \code{SuperSurv}.
#' @param newdata A data.frame of new covariate values.
#' @param new.times A numeric vector of times at which to predict survival.
#' @param onlySL Logical. If TRUE, only uses models with weights > threshold.
#' @param threshold Numeric. The weight threshold for onlySL.
#' @param ... Additional ignored arguments.
#' @return A list containing:
#' \itemize{
#'   \item \code{event.SL.predict}: A numeric matrix (rows = observations, columns = times) 
#'     of the final predicted survival probabilities from the ensemble.
#'   \item \code{event.library.predict}: A 3D numeric array (observations x times x models) 
#'     containing the individual survival predictions from each base learner.
#'   \item \code{cens.SL.predict}: A numeric matrix of the predicted censoring probabilities.
#'   \item \code{cens.library.predict}: A 3D numeric array of the individual censoring predictions.
#' }
#' @export
predict.SuperSurv <- function (object, newdata, new.times, onlySL = FALSE, threshold = 1e-4, ...) {
  
  # 1. Return training predictions if no new data is provided
  if (missing(newdata)) {
    out <- list(
      event.SL.predict      = object$event.SL.predict,
      cens.SL.predict       = object$cens.SL.predict,
      event.library.predict = object$event.library.predict,
      cens.library.predict  = object$cens.library.predict
    )
    return(out)
  }
  
  if (missing(new.times)) stop("new.times must be specified for new predictions.")
  
  # 2. Safety Check
  if (!object$control$saveFitLibrary) {
    stop("This SuperSurv fit was created using control$saveFitLibrary = FALSE. New predictions cannot be made.")
  }
  
  # 3. Setup Arrays
  event.k <- nrow(object$event.libraryNames)
  event.pred <- array(0, dim=c(nrow(newdata), length(new.times), event.k))
  dimnames(event.pred)[[3]] <- apply(object$event.libraryNames, 1, paste, collapse = '_')
  dimnames(event.pred)[[2]] <- new.times
  
  cens.k <- nrow(object$cens.libraryNames)
  cens.pred <- array(0, dim=c(nrow(newdata), length(new.times), cens.k))
  dimnames(cens.pred)[[3]] <- apply(object$cens.libraryNames, 1, paste, collapse = '_')
  dimnames(cens.pred)[[2]] <- new.times
  
  # 4. Filter by Threshold (onlySL)
  if (onlySL) {
    event.whichLibrary <- which(object$event.coef > threshold)
    if(length(event.whichLibrary) == 0) event.whichLibrary <- which(object$event.coef > 0) # Safety
    event.coef <- object$event.coef
    event.coef[-event.whichLibrary] <- 0
    event.coef <- event.coef / sum(event.coef)
    
    cens.whichLibrary <- which(object$cens.coef > threshold)
    if(length(cens.whichLibrary) == 0) cens.whichLibrary <- which(object$cens.coef > 0) # Safety
    cens.coef <- object$cens.coef
    cens.coef[-cens.whichLibrary] <- 0
    cens.coef <- cens.coef / sum(cens.coef)
  } else {
    event.whichLibrary <- seq(event.k)
    event.coef <- object$event.coef
    cens.whichLibrary <- seq(cens.k)
    cens.coef <- object$cens.coef
  }
  
  # ----------------------------------------------------------------------------
  # 5. Predict Event Models
  # ----------------------------------------------------------------------------
  for (mm in event.whichLibrary) {
    # Apply Variable Screening (if any)
    newdataMM <- subset(newdata, select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
    
    # Call our specific base learner prediction wrappers
    event.pred[, , mm] <- do.call("predict", list(
      object    = object$event.fitLibrary[[mm]],
      newX      = newdataMM, 
      new.times = new.times,
      ...
    ))
  }
  
  # Combine Event Learners (Using the 2D matrix force fix!)
  event.SL.predict <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(new.times))
  K_event <- length(event.coef)
  
  for (j in seq_along(new.times)) {
    tmp_mat <- matrix(event.pred[, j, , drop = FALSE], nrow = nrow(newdata), ncol = K_event)
    event.SL.predict[, j] <- tmp_mat %*% event.coef
  }
  
  # ----------------------------------------------------------------------------
  # 6. Predict Censoring Models
  # ----------------------------------------------------------------------------
  for (mm in cens.whichLibrary) {
    # Apply Variable Screening (if any)
    newdataMM <- subset(newdata, select = object$cens.whichScreen[object$cens.SL.library$library[mm, 2], ])
    
    # Call our specific base learner prediction wrappers
    cens.pred[, , mm] <- do.call("predict", list(
      object    = object$cens.fitLibrary[[mm]],
      newX      = newdataMM,
      new.times = new.times,
      ...
    ))
  }
  
  # Combine Censoring Learners
  cens.SL.predict <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(new.times))
  K_cens <- length(cens.coef)
  
  for (j in seq_along(new.times)) {
    tmp_mat <- matrix(cens.pred[, j, , drop = FALSE], nrow = nrow(newdata), ncol = K_cens)
    cens.SL.predict[, j] <- tmp_mat %*% cens.coef
  }
  
  # Ensure strict probability bounds to fix floating-point drift
  event.SL.predict <- pmax(pmin(event.SL.predict, 1), 0)
  cens.SL.predict <- pmax(pmin(cens.SL.predict, 1), 0)
  
  # ----------------------------------------------------------------------------
  # 7. Output
  # ----------------------------------------------------------------------------
  out <- list(
    event.SL.predict = event.SL.predict, 
    event.library.predict = event.pred,
    cens.SL.predict = cens.SL.predict, 
    cens.library.predict = cens.pred
  )
  return(out)
}