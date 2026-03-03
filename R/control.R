#' Control parameters for the SuperSurv Ensemble
#'
#' @param max.SL.iter Maximum iterations for the iterative weighting algorithm. Default 20.
#' @param event.t.grid Optional time grid for event risk calculation.
#' @param cens.t.grid Optional time grid for censoring risk calculation.
#' @param saveFitLibrary Logical. If TRUE (default), saves models for future predictions.
#' @param initWeightAlg The learner used for the very first step of IPCW.
#' @param initWeight Whether to start by fitting "censoring" or "event" weights.
#' @param ... Additional ignored arguments.
#' @return A list of control parameters.
#' @export
SuperSurv.control <- function(max.SL.iter = 20,
                              event.t.grid = NULL,
                              cens.t.grid = NULL,
                              saveFitLibrary = TRUE,
                              initWeightAlg = "surv.coxph", # Changed to a faster default
                              initWeight = "censoring",
                              ...) {

  if(!(initWeight %in% c("event", "censoring"))) {
    stop("initWeight must be one of 'event' or 'censoring'.")
  }

  # We return a list. If t.grids are NULL, the main SuperSurv function
  # will calculate them based on the observed data.
  list(
    max.SL.iter = max.SL.iter,
    event.t.grid = event.t.grid,
    cens.t.grid = cens.t.grid,
    saveFitLibrary = saveFitLibrary,
    initWeightAlg = initWeightAlg,
    initWeight = initWeight
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
#' @export
SuperSurv.CV.control <- function(V = 10L,
                                 stratifyCV = TRUE,
                                 shuffle = TRUE,
                                 validRows = NULL) {

  V <- as.integer(V)

  # Logic check: Can't have more folds than observations
  # (This will be double-checked in the main function once N is known)

  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("validRows must be a list of length V.")
    }
    if (!identical(V, length(validRows))) {
      # If user provided 5 folds of indices but set V=10, we must stop.
      V <- length(validRows)
      warning(paste("V adjusted to match length of validRows:", V))
    }
  }

  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}
