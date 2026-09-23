#' Control parameters for the SuperSurv Ensemble
#'
#' @param max.SL.iter Maximum iterations for the iterative weighting algorithm. Default 20.
#' @param event.t.grid Optional time grid for event risk calculation.
#' @param cens.t.grid Optional time grid for censoring risk calculation.
#' @param saveFitLibrary Logical. If TRUE (default), saves models for future predictions.
#' @param initWeightAlg The learner used for the very first step of IPCW.
#' @param initWeight Whether to start by fitting "censoring" or "event" weights.
#' @param tol Positive convergence tolerance for the maximum change in the
#'   out-of-fold event and censoring ensemble predictions.
#' @param traceIter Logical. If TRUE, reports iteration diagnostics.
#' @param ipcw.floor Smallest censoring or event survival probability used in
#'   an IPCW denominator.
#' @param ipcw.cap Largest inverse-probability weight used in an IPCW loss.
#' @param logloss.eps Probability clipping constant used only while evaluating
#'   logarithms in the IPCW log-loss.
#' @param optimizer.tol Positive convergence tolerance for metalearner
#'   optimization.
#' @param optimizer.maxit Maximum number of metalearner optimization
#'   iterations.
#' @param truncation.warn.fraction Fraction of IPCW rows stabilized by flooring
#'   or truncation above which a warning is issued.
#' @return A list of control parameters.
#' @keywords internal
SuperSurv.control <- function(max.SL.iter = 20,
                              event.t.grid = NULL,
                              cens.t.grid = NULL,
                              saveFitLibrary = TRUE,
                              initWeightAlg = "surv.coxph",
                              initWeight = "censoring",
                              tol = 1e-5,
                              traceIter = FALSE,
                              ipcw.floor = 1e-4,
                              ipcw.cap = 100,
                              logloss.eps = 1e-10,
                              optimizer.tol = 1e-8,
                              optimizer.maxit = 10000L,
                              truncation.warn.fraction = 0.05) {

  assert_scalar_numeric <- function(x, name, lower = -Inf, upper = Inf,
                                    lower_open = FALSE, upper_open = FALSE) {
    valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
    if (valid) {
      valid <- if (lower_open) x > lower else x >= lower
      valid <- valid && if (upper_open) x < upper else x <= upper
    }
    if (!valid) {
      interval <- paste0(
        if (lower_open) "(" else "[", lower, ", ", upper,
        if (upper_open) ")" else "]"
      )
      stop("`", name, "` must be a single finite number in ", interval, ".",
           call. = FALSE)
    }
  }

  if (!is.numeric(max.SL.iter) || length(max.SL.iter) != 1L ||
      is.na(max.SL.iter) || max.SL.iter < 1 || max.SL.iter != as.integer(max.SL.iter)) {
    stop("`max.SL.iter` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(optimizer.maxit) || length(optimizer.maxit) != 1L ||
      is.na(optimizer.maxit) || optimizer.maxit < 1 ||
      optimizer.maxit != as.integer(optimizer.maxit)) {
    stop("`optimizer.maxit` must be a positive integer.", call. = FALSE)
  }

  for (name in c("saveFitLibrary", "traceIter")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  if (!is.character(initWeightAlg) || length(initWeightAlg) != 1L ||
      is.na(initWeightAlg) || !nzchar(initWeightAlg)) {
    stop("`initWeightAlg` must be one non-empty learner function name.",
         call. = FALSE)
  }

  if (!is.character(initWeight) || length(initWeight) != 1L ||
      is.na(initWeight) || !(initWeight %in% c("event", "censoring"))) {
    stop("`initWeight` must be either 'event' or 'censoring'.", call. = FALSE)
  }

  assert_scalar_numeric(tol, "tol", lower = 0, lower_open = TRUE)
  assert_scalar_numeric(ipcw.floor, "ipcw.floor", lower = 0, upper = 1,
                        lower_open = TRUE, upper_open = TRUE)
  assert_scalar_numeric(ipcw.cap, "ipcw.cap", lower = 1)
  assert_scalar_numeric(logloss.eps, "logloss.eps", lower = 0, upper = 0.5,
                        lower_open = TRUE, upper_open = TRUE)
  assert_scalar_numeric(optimizer.tol, "optimizer.tol", lower = 0,
                        lower_open = TRUE)
  assert_scalar_numeric(truncation.warn.fraction, "truncation.warn.fraction",
                        lower = 0, upper = 1)

  list(
    max.SL.iter = as.integer(max.SL.iter),
    event.t.grid = event.t.grid,
    cens.t.grid = cens.t.grid,
    saveFitLibrary = saveFitLibrary,
    initWeightAlg = initWeightAlg,
    initWeight = initWeight,
    tol = tol,
    traceIter = traceIter,
    ipcw.floor = ipcw.floor,
    ipcw.cap = ipcw.cap,
    logloss.eps = logloss.eps,
    optimizer.tol = optimizer.tol,
    optimizer.maxit = as.integer(optimizer.maxit),
    truncation.warn.fraction = truncation.warn.fraction,
    time.weighting = "uniform"
  )
}





#' Control parameters for Cross-Validation in SuperSurv
#'
#' @param V Number of folds. Default is 10.
#' @param stratifyCV Logical. If TRUE, ensures event rates are balanced across folds.
#' @param shuffle Logical. If TRUE, shuffles rows before splitting.
#' @param validRows Optional custom list of indices for folds.
#'
#' @return A list of CV parameters.
#' @keywords internal
SuperSurv.CV.control <- function(V = 10L,
                                 stratifyCV = TRUE,
                                 shuffle = TRUE,
                                 validRows = NULL) {
  if (!is.numeric(V) || length(V) != 1L || !is.finite(V) ||
      V != as.integer(V) || V < 2L) {
    stop("`V` must be an integer of at least 2.", call. = FALSE)
  }
  V <- as.integer(V)
  .validate_scalar_logical(stratifyCV, "stratifyCV")
  .validate_scalar_logical(shuffle, "shuffle")

  # Logic check: Can't have more folds than observations
  # (This will be double-checked in the main function once N is known)

  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("`validRows` must be NULL or a list with one element per fold.",
           call. = FALSE)
    }
    if (!identical(V, length(validRows))) {
      # If user provided 5 folds of indices but set V=10, we must stop.
      V <- length(validRows)
      warning(paste("V adjusted to match length of validRows:", V))
    }
  }

  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}
