#' Predict method for SuperSurv fits
#'
#' Obtains predicted survival probabilities from a fitted SuperSurv ensemble.
#'
#' @param object A fitted object of class \code{SuperSurv}.
#' @param newdata A data.frame of new covariate values.
#' @param new.times A numeric vector of times at which to predict survival.
#' @param type Character string specifying the prediction output. Use
#'   \code{"event"} for the event survival matrix, \code{"censoring"} for the
#'   censoring survival matrix, or \code{"both"} for the full list of outputs.
#' @param onlySL Logical. If TRUE, only uses models with weights > threshold.
#' @param threshold Numeric. The weight threshold for onlySL.
#' @param ... Additional ignored arguments.
#' @return If \code{type = "event"} or \code{type = "censoring"}, a numeric
#'   matrix with rows corresponding to observations and columns corresponding to
#'   \code{new.times}. If \code{type = "both"}, a list containing:
#' \itemize{
#'   \item \code{event.predict}: A numeric matrix of final event survival predictions.
#'   \item \code{event.library.predict}: A 3D numeric array of event learner predictions.
#'   \item \code{cens.predict}: A numeric matrix of final censoring survival predictions.
#'   \item \code{cens.library.predict}: A 3D numeric array of censoring learner predictions.
#' }
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   newX <- X[1:10, , drop = FALSE]
#'   new.times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = new.times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   pred_event <- predict(
#'     object = fit,
#'     newdata = newX,
#'     new.times = new.times,
#'     type = "event"
#'   )
#'
#'   dim(pred_event)
#' }
#' @export
predict.SuperSurv <- function (object, newdata, new.times,
                               type = c("both", "event", "censoring"),
                               onlySL = FALSE, threshold = 1e-4, ...) {
  .validate_SuperSurv_object(object)
  type <- match.arg(type)
  .validate_scalar_logical(onlySL, "onlySL")
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("`threshold` must be one finite, non-negative number.", call. = FALSE)
  }

  # 1. Return training predictions if no new data is provided
  if (missing(newdata)) {
    if (!missing(new.times)) {
      stop("Supply `newdata` when requesting predictions at `new.times`.",
           call. = FALSE)
    }
    if (is.null(object$event.predict) || is.null(object$cens.predict)) {
      stop(
        "This SuperSurv object was fitted without stored prediction inputs. Supply both `newdata` and `new.times`.",
        call. = FALSE
      )
    }
    out <- list(
      event.predict      = object$event.predict,
      cens.predict       = object$cens.predict,
      event.library.predict = object$event.library.predict,
      cens.library.predict  = object$cens.library.predict
    )
    return(.predict_SuperSurv_output(out, type = type))
  }

  if (missing(new.times)) {
    stop("`new.times` must be specified for new predictions.", call. = FALSE)
  }
  new.times <- .validate_time_grid(new.times, "new.times")
  newdata <- .validate_prediction_newdata(newdata, object)

  # 2. Safety Check
  if (!isTRUE(object$control$saveFitLibrary)) {
    stop("This SuperSurv fit was created using `control$saveFitLibrary = FALSE`; refit with saved learners to predict on new data.",
         call. = FALSE)
  }
  if (!is.list(object$event.fitLibrary) || !is.list(object$cens.fitLibrary)) {
    stop("`object` does not contain valid saved event and censoring learner fits.",
         call. = FALSE)
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
    newdataMM <- subset(newdata, select = object$event.whichScreen[object$event.library$library[mm, 2], ])

    # Call our specific base learner prediction wrappers
    learner_name <- dimnames(event.pred)[[3L]][mm]
    learner_prediction <- tryCatch(
      do.call("predict", list(
        object = object$event.fitLibrary[[mm]],
        newdata = newdataMM,
        new.times = new.times,
        ...
      )),
      error = function(error) {
        stop("Prediction from event learner '", learner_name, "' failed: ",
             conditionMessage(error), call. = FALSE)
      }
    )
    event.pred[, , mm] <- .validate_learner_output(
      list(pred = learner_prediction),
      learner = learner_name,
      n_observations = nrow(newdata),
      times = new.times,
      context = "prediction on `newdata`"
    )
  }

  # Combine Event Learners (Using the 2D matrix force fix!)
  event.predict <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(new.times))
  K_event <- length(event.coef)

  for (j in seq_along(new.times)) {
    tmp_mat <- matrix(event.pred[, j, , drop = FALSE], nrow = nrow(newdata), ncol = K_event)
    event.predict[, j] <- tmp_mat %*% event.coef
  }

  if (type == "event") return(pmin(pmax(event.predict, 0), 1))

  # ----------------------------------------------------------------------------
  # 6. Predict Censoring Models
  # ----------------------------------------------------------------------------
  for (mm in cens.whichLibrary) {
    # Apply Variable Screening (if any)
    newdataMM <- subset(newdata, select = object$cens.whichScreen[object$cens.library$library[mm, 2], ])

    # Call our specific base learner prediction wrappers
    learner_name <- dimnames(cens.pred)[[3L]][mm]
    learner_prediction <- tryCatch(
      do.call("predict", list(
        object = object$cens.fitLibrary[[mm]],
        newdata = newdataMM,
        new.times = new.times,
        ...
      )),
      error = function(error) {
        stop("Prediction from censoring learner '", learner_name, "' failed: ",
             conditionMessage(error), call. = FALSE)
      }
    )
    cens.pred[, , mm] <- .validate_learner_output(
      list(pred = learner_prediction),
      learner = learner_name,
      n_observations = nrow(newdata),
      times = new.times,
      context = "prediction on `newdata`"
    )
  }

  # Combine Censoring Learners
  cens.predict <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(new.times))
  K_cens <- length(cens.coef)

  for (j in seq_along(new.times)) {
    tmp_mat <- matrix(cens.pred[, j, , drop = FALSE], nrow = nrow(newdata), ncol = K_cens)
    cens.predict[, j] <- tmp_mat %*% cens.coef
  }

  # Ensure strict probability bounds to fix floating-point drift
  event.predict <- pmax(pmin(event.predict, 1), 0)
  cens.predict <- pmax(pmin(cens.predict, 1), 0)

  # ----------------------------------------------------------------------------
  # 7. Output
  # ----------------------------------------------------------------------------
  out <- list(
    event.predict = event.predict,
    event.library.predict = event.pred,
    cens.predict = cens.predict,
    cens.library.predict = cens.pred
  )
  return(.predict_SuperSurv_output(out, type = type))
}

.predict_SuperSurv_output <- function(out, type) {
  switch(
    type,
    both = out,
    event = out$event.predict,
    censoring = out$cens.predict
  )
}
