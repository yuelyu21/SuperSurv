#' Super Learner for conditional survival functions
#'
#' Orchestrates the cross-validation, metalearner optimization, and prediction 
#' for an ensemble of survival base learners.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame for prediction (defaults to X).
#' @param new.times Times at which to obtain predicted survivals.
#' @param event.SL.library Character vector of prediction algorithms for the event.
#' @param cens.SL.library Character vector of prediction algorithms for censoring.
#' @param id Cluster identification variable.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param control List of control parameters for the Super Learner.
#' @param cvControl List of control parameters for cross-validation.
#' @param obsWeights Observation weights.
#' @param metalearner Character string specifying the optimizer (e.g., "least_squares").
#' @param parallel Logical. If TRUE, uses future.apply for parallel execution.
#' @param discrete Logical. If TRUE, selects the single best learner (weight = 1).
#'                 If FALSE (default), computes the optimal ensemble weights.
#' @param nFolds Number of cross-validation folds (default: 10).
#' @return A list of class \code{SuperSurv} containing:
#' \itemize{
#'   \item \code{call}: The matched function call.
#'   \item \code{event.SL.predict}: Matrix of in-sample cross-validated survival predictions.
#'   \item \code{cens.SL.predict}: Matrix of in-sample cross-validated censoring predictions.
#'   \item \code{event.coef}: Numeric vector of optimized ensemble weights for the event.
#'   \item \code{cens.coef}: Numeric vector of optimized ensemble weights for censoring.
#'   \item \code{event.library.predict}: 3D array of cross-validated predictions from individual event learners.
#'   \item \code{event.libraryNames}: Data frame detailing the algorithms and screeners used.
#'   \item \code{event.fitLibrary}: List of the fitted base learner models (if \code{saveFitLibrary = TRUE}).
#'   \item \code{times}: The time grid used for evaluation.
#' }
#' @export
SuperSurv <- function(time, event, X, newX = NULL, new.times, 
                      event.SL.library, cens.SL.library, 
                      id = NULL, verbose = FALSE, 
                      control = list(), cvControl = list(), 
                      obsWeights = NULL, metalearner = "least_squares",
                      parallel = FALSE, discrete = FALSE, nFolds = 10)  {
  
  # ----------------------------------------------------------------------------
  # 1. Setup & Checks
  # ----------------------------------------------------------------------------
  time_start <- proc.time()
  
  if (!is.null(dim(time)) && ncol(time) > 0) stop("time must be an (n x 1) numeric vector.")
  if (!is.null(dim(event)) && ncol(event) > 0) stop("event must be an (n x 1) numeric vector.")
  time <- as.numeric(time)
  event <- as.numeric(event)
  if (is.null(obsWeights)) obsWeights <- rep(1, length(time))
  if (is.null(newX)) newX <- X
  
  varNames <- colnames(X)
  N <- dim(X)[1L]
  p <- dim(X)[2L]
  
  # Internal Checker (from internals.R)
  .checkInputs(time=time, event=event, X=X, newX=newX, id=id, obsWeights=obsWeights, verbose=verbose)
  
  # ----------------------------------------------------------------------------
  # 2. Parallel Setup
  # ----------------------------------------------------------------------------
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required for parallel = TRUE. Please install it.")
    }
    lapply_fun <- future.apply::future_lapply
    parallel_args <- list(future.seed = TRUE) 
  } else {
    lapply_fun <- lapply
    parallel_args <- list()
  }
  
  # ----------------------------------------------------------------------------
  # 3. Time Grid Setup (Event) - With Memory Fix
  # ----------------------------------------------------------------------------
  if (is.null(control$event.t.grid)) {
    control$event.t.grid <- seq(0, max(time[event == 1]), length.out = 250)
  } else {
    if (!is.null(dim(control$event.t.grid)) && ncol(control$event.t.grid) > 0) 
      stop("control$event.t.grid must be a numeric vector.")
    control$event.t.grid <- sort(unique(as.numeric(control$event.t.grid)))
    if(any(is.na(control$event.t.grid))) stop("No missing values allowed in event.t.grid")
    if (any(control$event.t.grid < 0)) stop("Values in event.t.grid must be non-negative")
    if (any(control$event.t.grid > max(time))) stop("Values in event.t.grid must not exceed max(time)")
  }
  
  if(length(control$event.t.grid) > 500) {
    if(verbose) warning("event.t.grid > 500 points. Reducing to quantiles for memory safety.")
    control$event.t.grid <- quantile(control$event.t.grid, probs = seq(0, 1, length.out = 500))
  }
  
  # ----------------------------------------------------------------------------
  # 4. Time Grid Setup (Censoring) - With Left-Continuity & Memory Fix
  # ----------------------------------------------------------------------------
  if (is.null(control$cens.t.grid)) {
    control$cens.t.grid <- seq(0, max(time[event == 0]), length.out = 250)
    epsilon <- max(min(diff(sort(unique(time)))) / 2, 1e-5)
    control$cens.t.grid <- control$cens.t.grid - epsilon
    control$cens.t.grid <- c(0, control$cens.t.grid[control$cens.t.grid > 0])
  } else {
    control$cens.t.grid <- sort(unique(as.numeric(control$cens.t.grid)))
    if(any(is.na(control$cens.t.grid))) stop("No missing values allowed in cens.t.grid")
    if (any(control$cens.t.grid < 0)) stop("Values in cens.t.grid must be non-negative")
    if (any(control$cens.t.grid > max(time))) stop("Values in cens.t.grid must not exceed max(time)")
  }
  
  if(length(control$cens.t.grid) > 500) {
    if(verbose) warning("cens.t.grid > 500 points. Reducing to quantiles for memory safety.")
    control$cens.t.grid <- quantile(control$cens.t.grid, probs = seq(0, 1, length.out = 500))
  }
  
  # ----------------------------------------------------------------------------
  # 5. Libraries & Folds
  # ----------------------------------------------------------------------------
  if (missing(cvControl)) cvControl <- list() 
  cvControl$V <- nFolds
  
  control <- do.call("SuperSurv.control", control)
  cvControl <- do.call("SuperSurv.CV.control", cvControl)
  
  event.library <- .createLibrary(event.SL.library)
  cens.library <- .createLibrary(cens.SL.library)
  
  event.k <- nrow(event.library$library)
  cens.k <- nrow(cens.library$library)
  event.kScreen <- length(event.library$screenAlgorithm)
  cens.kScreen <- length(cens.library$screenAlgorithm)
  
  event.Z <- array(NA, dim = c(N, length(control$event.t.grid), event.k))
  cens.Z <- array(NA, dim = c(N, length(control$cens.t.grid), cens.k))
  
  event.libraryNames <- cbind(predAlgorithm = event.library$library$predAlgorithm,
                              screenAlgorithm = event.library$screenAlgorithm[event.library$library$rowScreen])
  cens.libraryNames <- cbind(predAlgorithm = cens.library$library$predAlgorithm,
                             screenAlgorithm = cens.library$screenAlgorithm[cens.library$library$rowScreen])
  
  if (p < 2 & !identical(event.library$screenAlgorithm, "All")) warning("Screening with p < 2.")
  if (p < 2 & !identical(cens.library$screenAlgorithm, "All")) warning("Screening with p < 2.")
  
  validRows <- .survCVFolds(N = N, id = id, event = event, cvControl = cvControl)
  if (is.null(id)) id <- seq(N)
  
  # ----------------------------------------------------------------------------
  # 6. Cross-Validation (Parallelized)
  # ----------------------------------------------------------------------------
  time_train_start <- proc.time()
  
  # --- Event CV ---
  if(verbose) message("Running Cross-Validation for Event Models...")
  
  event_args <- c(list(X = validRows, FUN = .crossValFUN, 
                       validRows = validRows,
                       time = time, event = event, dataX = X, id = id, obsWeights = obsWeights,
                       t.grid = control$event.t.grid, library = event.library,
                       kScreen = event.kScreen, k = event.k, p = p, verbose = verbose),
                  parallel_args) 
  
  event.crossValFUN_out <- do.call(lapply_fun, event_args)
  
  for(v in 1:cvControl$V) {
    event.Z[validRows[[v]], ,] <- event.crossValFUN_out[[v]]$out
  }
  
  event.errorsInCVLibrary <- apply(event.Z, 3, anyNA)
  if (sum(event.errorsInCVLibrary) > 0) event.Z[, , as.logical(event.errorsInCVLibrary)] <- 0
  if (all(event.Z == 0)) stop("All algorithms dropped from event library")
  
  # --- Censoring CV ---
  if(verbose) message("Running Cross-Validation for Censoring Models...")
  
  cens_args <- c(list(X = validRows, FUN = .crossValFUN,
                      validRows = validRows,
                      time = time, event = 1-event, dataX = X, id = id, obsWeights = obsWeights,
                      t.grid = control$cens.t.grid, library = cens.library,
                      kScreen = cens.kScreen, k = cens.k, p = p, verbose = verbose),
                 parallel_args)
  
  cens.crossValFUN_out <- do.call(lapply_fun, cens_args)
  
  for(v in 1:cvControl$V) {
    cens.Z[validRows[[v]], ,] <- cens.crossValFUN_out[[v]]$out
  }
  
  cens.errorsInCVLibrary <- apply(cens.Z, 3, anyNA)
  if (sum(cens.errorsInCVLibrary) > 0) cens.Z[, , as.logical(cens.errorsInCVLibrary)] <- 0
  if (all(cens.Z == 0)) stop("All algorithms dropped from censoring library")
  
  # ----------------------------------------------------------------------------
  # 7. Metalearner (With Selection)
  # ----------------------------------------------------------------------------
  getCoef <- .surviterativeSL(event.Z = event.Z, cens.Z = cens.Z, time = time, event = event, X = X,
                              obsWeights = obsWeights, id = id, control = control, verbose = verbose,
                              event.errorsInLibrary = event.errorsInCVLibrary, 
                              cens.errorsInLibrary = cens.errorsInCVLibrary,
                              metalearner = metalearner) 
  
  # if(verbose) {
  #   message("DEBUG: Metalearner coefficients returned.")
  # }
  
  event.coef <- getCoef$event.coef
  cens.coef <- getCoef$cens.coef
  event.cvRisks <- getCoef$event.cvRisks
  cens.cvRisks <- getCoef$cens.cvRisks
  
  names(event.coef) <- names(event.cvRisks) <- apply(event.libraryNames, 1, paste, collapse = "_")
  names(cens.coef) <- names(cens.cvRisks) <- apply(cens.libraryNames, 1, paste, collapse = "_")
  
  time_train <- proc.time() - time_train_start
  
# ----------------------------------------------------------------------------
  # 8. Full Fits (Prediction on newX)
  # ----------------------------------------------------------------------------
  if(verbose) message("Fitting full models on entire training set...")
  time_predict_start <- proc.time()
  
  # --- SCREENING ---
  event.whichScreen <- do.call(rbind, lapply(event.library$screenAlgorithm, FUN = .screenFun,
                                             list = list(time = time, event = event, id = id, X = X, obsWeights = obsWeights)))
  
  cens.whichScreen <- do.call(rbind, lapply(cens.library$screenAlgorithm, FUN = .screenFun,
                                            list = list(time = time, event = 1 - event, id = id, X = X, obsWeights = obsWeights)))
  
  # --- EVENT MODELS FIT ---
  event.pred <- lapply(seq(event.k), FUN = .predFun,
                       lib = event.library$library, time = time, event = event,
                       dataX = X, newX = newX, t.grid = new.times,
                       whichScreen = event.whichScreen, id = id,
                       obsWeights = obsWeights, verbose = verbose, control = control,
                       libraryNames = event.libraryNames)
  
  # Initialize 3D Array: [Patients x Times x Models]
  event.libraryPred <- array(NA, dim = c(nrow(newX), length(new.times), event.k))
  
  # [SAFETY FIX] Loop with check to prevent crash on failed models
  for (j in 1:event.k) {
    if (!is.null(event.pred[[j]]$out)) {
      event.libraryPred[, , j] <- event.pred[[j]]$out
    }
  }
  
  # Flag models that produced NAs (or were NULL/skipped)
  event.errorsInLibrary <- apply(event.libraryPred, 3, anyNA)
  
  # --- CENSORING MODELS FIT ---
  cens.pred <- lapply(seq(cens.k), FUN = .predFun,
                      lib = cens.library$library, time = time, event = 1 - event,
                      dataX = X, newX = newX, t.grid = new.times,
                      whichScreen = cens.whichScreen, id = id,
                      obsWeights = obsWeights, verbose = verbose, control = control,
                      libraryNames = cens.libraryNames)
  
  cens.libraryPred <- array(NA, dim = c(nrow(newX), length(new.times), cens.k))
  
  for (j in 1:cens.k) {
    if (!is.null(cens.pred[[j]]$out)) {
      cens.libraryPred[, , j] <- cens.pred[[j]]$out
    }
  }
  
  cens.errorsInLibrary <- apply(cens.libraryPred, 3, anyNA)

  # ----------------------------------------------------------------------------
  # 8b. Discrete Super Learner Logic (Override Weights)
  # ----------------------------------------------------------------------------
  if (discrete) {
    if(verbose) message("Discrete SL: Selecting the single best learner...")
    
    # Event Winner
    best_idx_event <- which.min(event.cvRisks)
    event.coef <- rep(0, length(event.cvRisks))
    event.coef[best_idx_event] <- 1
    if(verbose) message(paste0("  Event Winner: ", names(event.cvRisks)[best_idx_event]))
    
    # Censoring Winner
    best_idx_cens <- which.min(cens.cvRisks)
    cens.coef <- rep(0, length(cens.cvRisks))
    cens.coef[best_idx_cens] <- 1
    if(verbose) message(paste0("  Censoring Winner: ", names(cens.cvRisks)[best_idx_cens]))
  }
  
  # ----------------------------------------------------------------------------
  # 9. Final Ensemble & Output
  # ----------------------------------------------------------------------------
  
  # --- EVENT ENSEMBLE ---
  if (event.k == 1) {
    # If only one model exists, just take its prediction matrix
    event.SL.predict <- event.libraryPred[, , 1]
  } else {
    event.SL.predict <- matrix(NA_real_, nrow = nrow(newX), ncol = length(new.times))
    
    # Identify which models in the library actually finished without error
    valid_idx <- !event.errorsInLibrary
    if (sum(valid_idx) == 0) stop("All event models failed prediction.")
    
    # Extract coefficients for ONLY the valid models and re-normalize if necessary
    # (Though NNLS usually handles this, we be safe)
    current_event_coef <- event.coef[valid_idx]
    
    for (j in seq_along(new.times)) {
      # [CRITICAL FIX]: Force to 2D matrix to prevent "non-conformable" errors
      # This handles cases where sum(valid_idx) is 1 or many.
      tmp_mat <- matrix(event.libraryPred[, j, valid_idx, drop = FALSE], 
                        nrow = nrow(newX), 
                        ncol = sum(valid_idx))
      
      event.SL.predict[, j] <- tmp_mat %*% current_event_coef
    }
  }
  
  # --- CENSORING ENSEMBLE ---
  if (cens.k == 1) {
    cens.SL.predict <- cens.libraryPred[, , 1]
  } else {
    cens.SL.predict <- matrix(NA_real_, nrow = nrow(newX), ncol = length(new.times))
    
    valid_idx <- !cens.errorsInLibrary
    if (sum(valid_idx) == 0) stop("All censoring models failed prediction.")
    
    current_cens_coef <- cens.coef[valid_idx]
    
    for (j in seq_along(new.times)) {
      # [CRITICAL FIX]: Force to 2D matrix
      tmp_mat <- matrix(cens.libraryPred[, j, valid_idx, drop = FALSE], 
                        nrow = nrow(newX), 
                        ncol = sum(valid_idx))
      
      cens.SL.predict[, j] <- tmp_mat %*% current_cens_coef
    }
  }
  
  # --- TIMING & METADATA ---
  time_predict <- proc.time() - time_predict_start
  time_end <- proc.time()
  
  times <- list(everything = time_end - time_start, 
                train = time_train,
                predict = time_predict)
  
  # --- CONSTRUCT FIT LIBRARIES ---
  if (control$saveFitLibrary) {
    # Extract the actual model objects from the wrapper output lists
    event.fitLibrary <- lapply(event.pred, "[[", "model_out")
    names(event.fitLibrary) <- apply(event.libraryNames, 1, paste, collapse = "_")
    
    cens.fitLibrary <- lapply(cens.pred, "[[", "model_out")
    names(cens.fitLibrary) <- apply(cens.libraryNames, 1, paste, collapse = "_")
  } else {
    event.fitLibrary <- NULL
    cens.fitLibrary <- NULL
  }
  
  # ----------------------------------------------------------------------------
  # 10. Return List (Full Diagnostic Suite)
  # ----------------------------------------------------------------------------
  out <- list(
    call = match.call(),
    event.SL.predict = event.SL.predict, 
    cens.SL.predict = cens.SL.predict,
    event.coef = event.coef, 
    cens.coef = cens.coef,
    event.library.predict = event.libraryPred, 
    cens.library.predict = cens.libraryPred,
    event.libraryNames = event.libraryNames, 
    cens.libraryNames = cens.libraryNames,
    event.SL.library = event.library, 
    cens.SL.library = cens.library,
    event.cvRisks = event.cvRisks, 
    cens.cvRisks = cens.cvRisks,
    event.errorsInCVLibrary = event.errorsInCVLibrary, 
    cens.errorsInCVLibrary = cens.errorsInCVLibrary,
    event.errorsInLibrary = event.errorsInLibrary,     
    cens.errorsInLibrary = cens.errorsInLibrary,
    event.Z = event.Z, 
    cens.Z = cens.Z,
    varNames = varNames, 
    validRows = validRows, 
    event.whichScreen = event.whichScreen, 
    cens.whichScreen = cens.whichScreen,
    event.fitLibrary = event.fitLibrary, 
    cens.fitLibrary = cens.fitLibrary,
    control = control, 
    cvControl = cvControl,
    times = times
  )
  
  class(out) <- "SuperSurv"
  return(out)
}