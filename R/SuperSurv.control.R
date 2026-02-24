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
