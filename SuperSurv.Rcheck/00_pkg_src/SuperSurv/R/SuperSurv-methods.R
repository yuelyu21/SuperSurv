#' Extract SuperSurv ensemble coefficients
#'
#' Returns the optimized convex-combination weights from a fitted
#' \code{SuperSurv} object.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param type Character string specifying which weights to return. Use
#'   \code{"event"} for the event-time ensemble, \code{"censoring"} for the
#'   censoring ensemble, or \code{"both"} for both sets of weights.
#' @param ... Additional arguments ignored.
#'
#' @return For \code{type = "event"} or \code{type = "censoring"}, a named
#'   numeric vector of ensemble weights. For \code{type = "both"}, a list with
#'   elements \code{event} and \code{censoring}.
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
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
#'   coef(fit)
#'   coef(fit, type = "both")
#' }
#' @method coef SuperSurv
#' @export
coef.SuperSurv <- function(object, type = c("event", "censoring", "both"), ...) {
  type <- match.arg(type)

  switch(
    type,
    event = object$event.coef,
    censoring = object$cens.coef,
    both = list(
      event = object$event.coef,
      censoring = object$cens.coef
    )
  )
}

#' Access SuperSurv ensemble weights
#'
#' Extracts the fitted event or censoring ensemble weights from a
#' \code{SuperSurv} object.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param ... Additional arguments ignored.
#'
#' @return A named numeric vector of ensemble weights.
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
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
#'   event_weights(fit)
#'   censor_weights(fit)
#' }
#' @export
event_weights <- function(object, ...) {
  UseMethod("event_weights")
}

#' @rdname event_weights
#' @method event_weights SuperSurv
#' @export
event_weights.SuperSurv <- function(object, ...) {
  coef(object, type = "event")
}

#' @rdname event_weights
#' @export
censor_weights <- function(object, ...) {
  UseMethod("censor_weights")
}

#' @rdname event_weights
#' @method censor_weights SuperSurv
#' @export
censor_weights.SuperSurv <- function(object, ...) {
  coef(object, type = "censoring")
}

#' Access SuperSurv learner names
#'
#' Returns the fitted learner names from a \code{SuperSurv} object.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param type Character string specifying whether to return event learner
#'   names, censoring learner names, or both.
#' @param ... Additional arguments ignored.
#'
#' @return For \code{type = "event"} or \code{type = "censoring"}, a character
#'   vector of learner names. For \code{type = "both"}, a list with elements
#'   \code{event} and \code{censoring}.
#'
#' @export
learner_names <- function(object, ...) {
  UseMethod("learner_names")
}

#' @rdname learner_names
#' @method learner_names SuperSurv
#' @export
learner_names.SuperSurv <- function(object,
                                    type = c("event", "censoring", "both"),
                                    ...) {
  type <- match.arg(type)

  event_names <- names(event_weights(object))
  censor_names <- names(censor_weights(object))

  switch(
    type,
    event = event_names,
    censoring = censor_names,
    both = list(event = event_names, censoring = censor_names)
  )
}

#' Access SuperSurv prediction evaluation times
#'
#' Returns the time grid used for the fitted \code{SuperSurv} predictions.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param ... Additional arguments ignored.
#'
#' @return A numeric vector of prediction evaluation times.
#'
#' @export
eval_times <- function(object, ...) {
  UseMethod("eval_times")
}

#' @rdname eval_times
#' @method eval_times SuperSurv
#' @export
eval_times.SuperSurv <- function(object, ...) {
  times <- object$eval.times
  if (is.null(times)) {
    pred_dimnames <- dimnames(object$event.predict)
    times <- suppressWarnings(as.numeric(pred_dimnames[[2L]]))
  }
  if (is.null(times) || length(times) == 0L || anyNA(times)) {
    stop("Prediction evaluation times are not available in this SuperSurv object.",
         call. = FALSE)
  }

  times
}

#' Access SuperSurv training variable names
#'
#' Returns the covariate names used to fit a \code{SuperSurv} object.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param ... Additional arguments ignored.
#'
#' @return A character vector of training variable names.
#'
#' @export
training_variables <- function(object, ...) {
  UseMethod("training_variables")
}

#' @rdname training_variables
#' @method training_variables SuperSurv
#' @export
training_variables.SuperSurv <- function(object, ...) {
  vars <- object$varNames
  if (is.null(vars)) {
    stop("Training variable names are not available in this SuperSurv object.", call. = FALSE)
  }

  vars
}

#' Access variables selected by SuperSurv screeners
#'
#' Returns the variables retained by a screening step for one or more fitted
#' event or censoring learners.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param type Character string specifying whether to inspect event or censoring
#'   learners.
#' @param learner Optional learner index or learner name. If omitted, selected
#'   variables are returned for every learner of the requested type.
#' @param ... Additional arguments ignored.
#'
#' @return A named list of character vectors, or a single character vector when
#'   \code{learner} has length one.
#'
#' @export
selected_variables <- function(object, ...) {
  UseMethod("selected_variables")
}

#' @rdname selected_variables
#' @method selected_variables SuperSurv
#' @export
selected_variables.SuperSurv <- function(object,
                                         type = c("event", "censoring"),
                                         learner = NULL,
                                         ...) {
  type <- match.arg(type)
  vars <- training_variables(object)

  screen <- switch(
    type,
    event = object$event.whichScreen,
    censoring = object$cens.whichScreen
  )

  names_by_learner <- learner_names(object, type = type)
  screen <- .SuperSurv_screen_by_learner(object, screen, type)

  if (is.null(learner)) {
    learner_idx <- seq_along(names_by_learner)
  } else if (is.character(learner)) {
    learner_idx <- match(learner, names_by_learner)
    if (anyNA(learner_idx)) {
      stop("Unknown learner name: ",
           paste(learner[is.na(learner_idx)], collapse = ", "),
           call. = FALSE)
    }
  } else {
    learner_idx <- as.integer(learner)
    if (anyNA(learner_idx) || any(learner_idx < 1L) ||
        any(learner_idx > length(names_by_learner))) {
      stop("'learner' must be a valid learner index or name.", call. = FALSE)
    }
  }

  out <- lapply(learner_idx, function(i) vars[as.logical(screen[i, ])])
  names(out) <- names_by_learner[learner_idx]

  if (length(out) == 1L) out[[1L]] else out
}

#' Print a SuperSurv fit
#'
#' Prints a concise description of a fitted \code{SuperSurv} object.
#'
#' @param x A fitted object of class \code{"SuperSurv"}.
#' @param digits Number of significant digits to use for displayed weights.
#' @param ... Additional arguments ignored.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @method print SuperSurv
#' @export
print.SuperSurv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("SuperSurv fit\n")
  cat("  Selection:", .SuperSurv_selection(x), "\n")
  cat("  Event learners:", .SuperSurv_n_learners(x$event.libraryNames), "\n")
  cat("  Censoring learners:", .SuperSurv_n_learners(x$cens.libraryNames), "\n")

  pred_dim <- dim(x$event.predict)
  if (length(pred_dim) == 2L) {
    cat("  Predictions:", pred_dim[1L], "observations x", pred_dim[2L], "times\n")
  }
  cat("  Evaluation times:", .SuperSurv_range_label(eval_times(x)), "\n")

  event_weights <- event_weights(x)
  active <- event_weights[event_weights > 0]
  if (length(active) > 0L) {
    active <- sort(active, decreasing = TRUE)
    cat("  Nonzero event weights:\n")
    print(round(active, digits))
  } else {
    cat("  Nonzero event weights: none\n")
  }

  invisible(x)
}

#' Summarize a SuperSurv fit
#'
#' Summarizes the fitted event and censoring ensembles, including learner names,
#' ensemble weights, cross-validated risks, error flags, prediction dimensions,
#' and recorded timing information.
#'
#' @param object A fitted object of class \code{"SuperSurv"}.
#' @param x A summary object produced by \code{summary.SuperSurv()}.
#' @param digits Number of significant digits to use for displayed weights and
#'   risks.
#' @param ... Additional arguments ignored.
#'
#' @return An object of class \code{"summary.SuperSurv"}, a list containing the
#'   matched call, selection mode, event and censoring learner summaries,
#'   prediction dimensions, and timing information.
#'
#' @method summary SuperSurv
#' @export
summary.SuperSurv <- function(object, ...) {
  out <- list(
    call = object$call,
    selection = .SuperSurv_selection(object),
    event = .SuperSurv_learner_summary(
      learner_names = object$event.libraryNames,
      weights = event_weights(object),
      cv_risks = object$event.cvRisks,
      cv_errors = object$event.errorsInCVLibrary,
      fit_errors = object$event.errorsInLibrary
    ),
    censoring = .SuperSurv_learner_summary(
      learner_names = object$cens.libraryNames,
      weights = censor_weights(object),
      cv_risks = object$cens.cvRisks,
      cv_errors = object$cens.errorsInCVLibrary,
      fit_errors = object$cens.errorsInLibrary
    ),
    prediction_dimensions = list(
      event = dim(object$event.predict),
      censoring = dim(object$cens.predict),
      event_library = dim(object$event.library.predict),
      censoring_library = dim(object$cens.library.predict)
    ),
    eval_times = eval_times(object),
    times = object$times
  )

  class(out) <- "summary.SuperSurv"
  out
}

#' @rdname summary.SuperSurv
#' @method print summary.SuperSurv
#' @export
print.summary.SuperSurv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Summary of SuperSurv fit\n")
  cat("  Selection:", x$selection, "\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Event ensemble:\n")
  print(.SuperSurv_format_summary_table(x$event, digits = digits), row.names = FALSE)
  cat("\nCensoring ensemble:\n")
  print(.SuperSurv_format_summary_table(x$censoring, digits = digits), row.names = FALSE)

  event_dim <- x$prediction_dimensions$event
  if (length(event_dim) == 2L) {
    cat("\nPredictions:", event_dim[1L], "observations x", event_dim[2L], "times\n")
  }
  if (!is.null(x$eval_times)) {
    cat("Evaluation times:", .SuperSurv_range_label(x$eval_times), "\n")
  }

  if (!is.null(x$times)) {
    elapsed <- vapply(x$times, function(z) unname(z[["elapsed"]]), numeric(1L))
    cat("Elapsed time (seconds):\n")
    print(round(elapsed, digits))
  }

  invisible(x)
}

.SuperSurv_selection <- function(object) {
  if (!is.null(object$call) && !is.null(object$call$selection)) {
    selection <- object$call$selection
    if (is.character(selection)) return(selection)
    return(paste(deparse(selection), collapse = " "))
  }

  "ensemble"
}

.SuperSurv_n_learners <- function(learner_names) {
  n <- nrow(learner_names)
  if (is.null(n)) length(learner_names) else n
}

.SuperSurv_learner_summary <- function(learner_names, weights, cv_risks,
                                       cv_errors, fit_errors) {
  learner_names <- as.data.frame(learner_names, stringsAsFactors = FALSE)
  if (!all(c("predAlgorithm", "screenAlgorithm") %in% names(learner_names))) {
    names(learner_names)[seq_len(min(2L, ncol(learner_names)))] <-
      c("predAlgorithm", "screenAlgorithm")[seq_len(min(2L, ncol(learner_names)))]
  }

  out <- data.frame(
    learner = names(weights),
    predAlgorithm = learner_names$predAlgorithm,
    screenAlgorithm = learner_names$screenAlgorithm,
    weight = unname(weights),
    risk = unname(cv_risks),
    status = .SuperSurv_learner_status(cv_errors, fit_errors),
    cv_error = as.logical(cv_errors),
    fit_error = as.logical(fit_errors),
    stringsAsFactors = FALSE
  )

  out[order(out$risk, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
}

.SuperSurv_format_summary_table <- function(x, digits) {
  x$weight <- round(x$weight, digits)
  x$risk <- round(x$risk, digits)
  x[, c("learner", "weight", "risk", "status"), drop = FALSE]
}

.SuperSurv_range_label <- function(x) {
  if (length(x) == 1L) return(format(x))
  paste0(length(x), " values from ", format(min(x)), " to ", format(max(x)))
}

.SuperSurv_learner_status <- function(cv_errors, fit_errors) {
  status <- rep("ok", length(cv_errors))
  status[as.logical(cv_errors)] <- "CV error"
  status[as.logical(fit_errors)] <- "fit error"
  status[as.logical(cv_errors) & as.logical(fit_errors)] <- "CV and fit error"
  status
}

.SuperSurv_screen_by_learner <- function(object, screen, type) {
  if (is.null(dim(screen))) {
    screen <- matrix(as.logical(screen), nrow = 1L)
  }

  n_learners <- length(learner_names(object, type = type))
  if (nrow(screen) == n_learners) {
    return(screen)
  }

  library_rows <- switch(
    type,
    event = object$event.library$library$rowScreen,
    censoring = object$cens.library$library$rowScreen
  )

  screen[library_rows, , drop = FALSE]
}
