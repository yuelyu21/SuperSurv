########## function previously inside the main function
### Get cross-validated survivals for S and T
#' @noRd
.crossValFUN <- function(valid, validRows, time, event, dataX, id, obsWeights, t.grid,
                         library, kScreen, k, p, verbose) {

  foldNum <- as.numeric(which(unlist(lapply(validRows, function(v) all.equal(v, valid) == TRUE))))
  tempLearn <- dataX[-valid, , drop = FALSE]
  tempTime <- time[-valid]
  tempEvent <- event[-valid]
  tempValid <- dataX[valid, , drop = FALSE]
  tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
  tempId <- id[-valid]
  tempObsWeights <- obsWeights[-valid]
  for (s in seq(kScreen)) {
    if(verbose) message(paste("CV ", library$screenAlgorithm[s],
                              ", fold ", foldNum, sep = ""))
    screen_fn <- get(library$screenAlgorithm[s])
    testScreen <- try(do.call(screen_fn, list(time = tempTime,
                                              event = tempEvent,
                                              X = tempLearn, id = tempId,
                                              obsWeights = tempObsWeights)))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,",
                    library$screenAlgorithm[s], ", with All()",
                    "\n "))
      tempWhichScreen[s, ] <- TRUE
    }
    else {
      tempWhichScreen[s, ] <- testScreen
    }
    if (verbose) {
      message(paste("Number of covariates in ", library$screenAlgorithm[s],
                    " is: ", sum(tempWhichScreen[s, ]), sep = ""))
    }
  }

  uniqueScreen <- unique(tempWhichScreen)
  screenMap <- apply(uniqueScreen, 1, function(row) which(apply(tempWhichScreen, 1, function(row2) all.equal(row, row2) == TRUE)))

  out <- array(NA, dim = c(nrow(tempValid), length(t.grid), k))

  for (predAlg in unique(library$library$predAlgorithm)) {
    if (verbose) message(paste("CV ", predAlg, ", fold ", foldNum, sep = ""))
    pred_fn <- get(predAlg)
    for(j in seq(nrow(uniqueScreen))) {
      testAlg <- try(do.call(pred_fn, list(time = tempTime, event = tempEvent,
                                           X = subset(tempLearn, select = uniqueScreen[j,], drop = FALSE),
                                           newX = subset(tempValid, select = uniqueScreen[j,], drop = FALSE),
                                           new.times = t.grid,
                                           id = tempId,
                                           obsWeights = tempObsWeights)))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", predAlg,
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      } else {
        libraryRows <- which(library$library$predAlgorithm == predAlg & library$library$rowScreen %in% unlist(screenMap[j]))
        for (row in libraryRows) {
          out[,,row] <- testAlg$pred
        }
      }
    }
  }


  invisible(list(out = out))
}




# Generate Full Screening logic
#' @noRd
.screenFun <- function(fun, list) {
  screen_fn <- get(fun)
  # Robust screen call
  testScreen <- try(do.call(screen_fn, list))
  if (inherits(testScreen, "try-error")) {
    warning(paste("Screening failed:", fun, "- defaulting to ALL"))
    out <- rep(TRUE, ncol(list$X))
  } else {
    out <- testScreen
  }
  return(out)
}


.predFun <- function(index, lib, time, event, dataX, newX, whichScreen, t.grid,
                     family, id, obsWeights, verbose, control, libraryNames) {
  if (verbose) {
    message(paste("full", libraryNames[index]))
  }
  pred_fn <- get(lib$predAlgorithm[index])
  testAlg <- try(do.call(pred_fn, list(time = time, event = event,
                                       X = subset(dataX, select = whichScreen[lib$rowScreen[index], ],
                                                  drop = FALSE),
                                       newX = subset(newX, select = whichScreen[lib$rowScreen[index],
                                       ], drop = FALSE), id = id,
                                       obsWeights = obsWeights, new.times = t.grid)))
  if (inherits(testAlg, "try-error")) {
    warning(paste("Error in algorithm", lib$predAlgorithm[index],
                  " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
    out <- rep.int(NA, times = nrow(newX))
    model_out <- NULL
  }
  else {
    out <- testAlg$pred
    if(control$saveFitLibrary) model_out <- testAlg$fit
    else model_out <- NULL
  }

  invisible(list(out = out, model_out = model_out))
}


########## function previously inside the main function
#' @noRd
.checkInputs <- function(time, event, X, newX, id, obsWeights, verbose) {
  if(any(time < 0)) stop("Only non-negative event/censoring times allowed!")
  if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
  if(any(is.na(time)) | any(is.na(event))) stop("No missing values allowed in time or event.")
  if(any(is.na(X)) | any(is.na(newX))) stop("No missing values allowed in X or new X.")
  if(length(time) != length(event) | length(time) != nrow(X)) stop("time and event must be n x 1 vectors and X must have n rows.")
  if(!is.data.frame(X) | !is.data.frame(newX)) stop("X and newX must be data frames.")
  if(!identical(names(X), names(newX))) stop("X and newX must have the same features.")
  if(any(obsWeights < 0)) stop("obsWeights < 0 not allowed.")
  if(!(verbose %in% c(TRUE, FALSE))) stop("verbose must be TRUE/FALSE.")
  if(!is.null(id) && !identical(length(id), length(time))) stop("id vector must have the same dimension as time")
}


#' @noRd
.createLibrary <- function (survSL.library)  {
  if (is.character(survSL.library)) {
    k <- length(survSL.library)
    whichScreen <- matrix(1, nrow = 1, ncol = k)
    screenAlgorithm <- "screen.all"
    library <- data.frame(predAlgorithm = survSL.library, rowScreen = 1,
                          stringsAsFactors = FALSE)
  }
  else if (is.list(survSL.library)) {
    predNames <- sapply(survSL.library, FUN = "[", 1)
    NumberScreen <- (sapply(survSL.library, FUN = length) - 1)
    if (sum(NumberScreen == 0) > 0) {
      for (ii in which(NumberScreen == 0)) {
        SL.library[[ii]] <- c(SL.library[[ii]], "screen.all")
        NumberScreen[ii] <- 1
      }
    }
    screenAlgorithmFull <- unlist(lapply(survSL.library, FUN = "[", -1))
    screenAlgorithm <- unique(screenAlgorithmFull)
    library <- data.frame(predAlgorithm = rep(predNames,
                                              times = NumberScreen), rowScreen = match(screenAlgorithmFull,
                                                                                       screenAlgorithm), stringsAsFactors = FALSE)
  }
  else {
    stop("format for survSL.library is not recognized")
  }
  out <- list(library = library, screenAlgorithm = screenAlgorithm)
  return(out)
}


#' @noRd
.survCVFolds <- function (N, id, event, cvControl) {
  if (!is.null(cvControl$validRows)) return(cvControl$validRows)
  stratifyCV <- cvControl$stratifyCV
  shuffle <- cvControl$shuffle
  V <- cvControl$V
  if (!stratifyCV) {
    if (shuffle) {
      if (is.null(id)) {
        validRows <- split(sample(1:N), rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(sample(1:n.id), rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
    else {
      if (is.null(id)) {
        validRows <- split(1:N, rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(1:n.id, rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
  }
  else {
    # if (sum(event) < V | sum(1-event) < V) {
    #   stop("number of (event = 1) or (event = 0) is less than the number of folds")
    # }
    if (shuffle) {
      if (is.null(id)) {
        event.0 <- which(event == 0)
        event.1 <- which(event == 1)
        rows.0 <- split(sample(event.0), rep(1:V, length = length(event.0)))
        rows.1 <- split(sample(event.1), rep(1:V, length = length(event.1)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          if (length(rows.0) >= vv) {
            if (length(rows.1) >= vv) validRows[[vv]] <- c(rows.0[[vv]], rows.1[[vv]])
            else validRows[[vv]] <- rows.0[[vv]]
          } else {
            validRows[[vv]] <- rows.1[[vv]]
          }
        }
      }
      else {
        stop("Stratified sampling with id not currently implemented. Either remove id or set control(stratifyCV = FALSE).")
      }
    }
    else {
      if (is.null(id)) {
        within.split <- suppressWarnings(tapply(1:N,
                                                INDEX = event, FUN = split, rep(1:V)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          validRows[[vv]] <- c(within.split[[1]][[vv]],
                               within.split[[2]][[vv]])
        }
      }
      else {
        stop("Stratified sampling with id not currently implemented. Either remove id or set control(stratifyCV = FALSE).")
      }
    }
  }
  return(validRows)
}



# ==============================================================================
# INTERNAL HELPER FUNCTIONS FOR SUPERSURV
# ==============================================================================

#' Internal function to perform iterative Super Learner optimization
#' @noRd
.surviterativeSL <- function(event.Z, cens.Z, time, event, X, obsWeights, id, control, verbose,
                             event.errorsInLibrary, cens.errorsInLibrary,
                             metalearner = "least_squares") {

  if(verbose) message("Performing iterative SuperLearner optimization...")

  event.k <- dim(event.Z)[3]
  cens.k <- dim(cens.Z)[3]
  N <- length(time)
  event.n.time <- length(control$event.t.grid)
  cens.n.time <- length(control$cens.t.grid)

  # Epsilon for offsets (Left-continuity)
  epsilon <- max(min(diff(sort(unique(time)))), 1e-5)

  # Flatten Arrays (Long Format)
  event.Z.long <- matrix(NA, nrow = N * event.n.time, ncol = event.k)
  for(j in seq(dim(event.Z)[3])) event.Z.long[,j] <- c(event.Z[,,j])

  cens.Z.long <- matrix(NA, nrow = N * cens.n.time, ncol = cens.k)
  for(j in seq(dim(cens.Z)[3])) cens.Z.long[,j] <- c(cens.Z[,,j])

  # Create "Observed" Predictions (Interpolated)
  event.Z.obs <- matrix(NA, nrow = N, ncol = event.k)
  for(i in seq(N)) {
    for(j in seq(event.k)) {
      event.Z.obs[i,j] <- stats::approx(control$event.t.grid, event.Z[i,,j], xout = time[i], method = 'constant', rule = 2, ties = mean)$y
    }
  }

  cens.Z.obs <- matrix(NA, nrow = N, ncol = cens.k)
  for(i in seq(N)) {
    for(j in seq(cens.k)) {
      cens.Z.obs[i,j] <- stats::approx(c(-1,control$cens.t.grid), c(1,cens.Z[i,,j]), xout = time[i] - epsilon, method = 'constant', rule = 2, ties = mean)$y
    }
  }

  # Setup Long Vectors for Optimization
  obsWeights.event.long <- rep(obsWeights, event.n.time)
  obsWeights.cens.long <- rep(obsWeights, cens.n.time)
  time.event.long <- rep(time, event.n.time)
  time.cens.long <- rep(time, cens.n.time)
  event.event.long <- rep(event, event.n.time)
  event.cens.long <- rep(event, cens.n.time)
  event.t.grid.long <- rep(control$event.t.grid, each = N)
  cens.t.grid.long <- rep(control$cens.t.grid, each = N)

  # Initial Weights
  initWeightAlg <- get(control$initWeightAlg)

  # Initialize variables
  obs.cens.vals <- NULL
  obs.event.vals <- NULL
  S.coef <- rep(0, event.k)
  G.coef <- rep(0, cens.k)

  # Initialization Step
  if (control$initWeight == "censoring") {
    initFit <- initWeightAlg(time = time, event = 1 - event, X = X, newX = X,
                             new.times = time - epsilon,
                             obsWeights = obsWeights, id = id)
    obs.cens.vals <- rep(diag(initFit$pred), length(control$event.t.grid))

    # Calculate Event Coefs
    S.coef[!event.errorsInLibrary] <- .survcomputeCoef(
      time = time.event.long, event = event.event.long,
      t.vals = event.t.grid.long, cens.vals = obs.cens.vals,
      preds = event.Z.long[,!event.errorsInLibrary, drop=FALSE],
      obsWeights = obsWeights.event.long,
      method = metalearner
    )
    obs.event.vals <- rep(c(event.Z.obs %*% S.coef), length(control$cens.t.grid))
  } else {
    initFit <- initWeightAlg(time = time, event = event, X = X, newX = X,
                             new.times = time,
                             obsWeights = obsWeights, id = id)
    obs.event.vals <- rep(diag(initFit$pred), length(control$cens.t.grid))
  }

  obs.event.vals[obs.event.vals == 0] <- min(obs.event.vals[obs.event.vals > 0])

  # Iteration Loop
  iter <- 1
  while(TRUE) {
    if(iter > control$max.SL.iter) {
      warning("Did not converge in ", control$max.SL.iter, " iterations")
      break
    }
    if(!is.null(obs.cens.vals)) obs.cens.vals.old <- obs.cens.vals
    if(!is.null(obs.event.vals)) obs.event.vals.old <- obs.event.vals

    # Update Censoring Weights
    G.coef[!cens.errorsInLibrary] <- .survcomputeCoef(
      time = time.cens.long, event = 1 - event.cens.long,
      t.vals = cens.t.grid.long, cens.vals = obs.event.vals,
      preds = cens.Z.long[,!cens.errorsInLibrary, drop=FALSE],
      obsWeights = obsWeights.cens.long,
      method = metalearner
    )
    obs.cens.vals <- rep(c(cens.Z.obs %*% G.coef), length(control$event.t.grid))
    obs.cens.vals[obs.cens.vals == 0] <- min(obs.cens.vals[obs.cens.vals > 0])

    # Update Event Weights
    S.coef[!event.errorsInLibrary] <- .survcomputeCoef(
      time = time.event.long, event = event.event.long,
      t.vals = event.t.grid.long, cens.vals = obs.cens.vals,
      preds = event.Z.long[,!event.errorsInLibrary, drop=FALSE],
      obsWeights = obsWeights.event.long,
      method = metalearner
    )
    obs.event.vals <- rep(c(event.Z.obs %*% S.coef), length(control$cens.t.grid))
    obs.event.vals[obs.event.vals == 0] <- min(obs.event.vals[obs.event.vals > 0])

    # Convergence Check
    if(!is.null(obs.cens.vals.old) & !is.null(obs.event.vals.old)) {
      cens.delta <- max(abs(obs.cens.vals - obs.cens.vals.old))
      event.delta <- max(abs(obs.event.vals - obs.event.vals.old))
      if(cens.delta + event.delta < 1e-5) {
        if(verbose) message("Converged in ", iter, " iterations.")
        break
      }
    }
    iter <- iter + 1
  }

  # ----------------------------------------------------------------------------
  # Calculate Final CV Risks (Dynamic based on Method)
  # ----------------------------------------------------------------------------
  calc_risk <- function(preds, time, event, t.grid, cens.vals, weights, method) {
    # Raw IPCW Pseudo-value
    Y_ipcw <- 1 - (as.numeric(time <= t.grid) * event / cens.vals)

    if (method == "least_squares") {
      # MSE can handle raw pseudo-values perfectly
      risks <- apply(preds, 2, function(col) {
        mean(weights * (Y_ipcw - col)^2, na.rm = TRUE)
      })
    } else {
      # Log-Loss REQUIRES the target to be bounded [0, 1]
      Y_clamp <- pmax(pmin(Y_ipcw, 1), 0)

      risks <- apply(preds, 2, function(col) {
        p <- pmax(pmin(col, 1 - 1e-15), 1e-15)
        loss <- -(weights * (Y_clamp * log(p) + (1 - Y_clamp) * log(1 - p)))
        mean(loss, na.rm = TRUE)
      })
    }
    return(risks)
  }


  event.cvRisks <- calc_risk(
    preds = event.Z.long, time = time.event.long, event = event.event.long,
    t.grid = event.t.grid.long, cens.vals = obs.cens.vals,
    weights = obsWeights.event.long, method = metalearner
  )


  cens.cvRisks <- calc_risk(
    preds = cens.Z.long, time = time.cens.long, event = event.cens.long,
    t.grid = cens.t.grid.long, cens.vals = obs.event.vals,
    weights = obsWeights.cens.long, method = metalearner
  )

  return(list(event.coef = S.coef, cens.coef = G.coef, event.cvRisks = event.cvRisks, cens.cvRisks = cens.cvRisks))
}

#' Internal function to compute ensemble weights
#' @importFrom nnls nnls
#' @noRd
.survcomputeCoef <- function(time, event, t.vals, cens.vals, preds, obsWeights, method = "least_squares") {

  # Safety: Remove NA columns
  valid_cols <- !apply(preds, 2, anyNA)
  if(sum(valid_cols) == 0) return(rep(1/ncol(preds), ncol(preds)))

  preds_clean <- preds[, valid_cols, drop=FALSE]
  k_clean <- ncol(preds_clean)
  if(k_clean == 1) {
    full_w <- rep(0, ncol(preds))
    full_w[valid_cols] <- 1
    return(full_w)
  }

  cens.vals[cens.vals < 1e-4] <- 1e-4

  # ----------------------------------------------------------------------------
  # METHOD 1: Least Squares (Brier Score)
  # ----------------------------------------------------------------------------
  if (method == "least_squares") {
    # IPCW "Observed" Outcome (1 if alive, 0 if dead, weighted)
    out <- 1 - as.numeric(time <= t.vals) * event / cens.vals

    fit.nnls <- nnls::nnls(sqrt(obsWeights) * preds_clean, sqrt(obsWeights) * out)
    coef <- coef(fit.nnls)
    if(sum(coef) == 0) coef <- rep(1, length(coef))

    final_w <- coef / sum(coef)
  }

  # ----------------------------------------------------------------------------
  # METHOD 2: Negative Log Likelihood (Log-Loss)
  # ----------------------------------------------------------------------------
  else if (method == "nloglik") {

    # Raw target
    Y_ipcw <- 1 - (as.numeric(time <= t.vals) * event / cens.vals)

    # CLAMP the target so cross-entropy doesn't diverge to negative infinity
    Y_clamp <- pmax(pmin(Y_ipcw, 1), 0)

    # Objective Function: Weighted Binary Cross Entropy
    fn <- function(w) {
      P_ens <- preds_clean %*% w
      P_ens <- pmax(pmin(P_ens, 1 - 1e-10), 1e-10)

      # Using the clamped target
      loss <- -(obsWeights * (Y_clamp * log(P_ens) + (1 - Y_clamp) * log(1 - P_ens)))
      return(mean(loss))
    }

    # Strategy: Optimize k weights, penalty for sum != 1
    fn_penalized <- function(w) {
      pen <- 1000 * (sum(w) - 1)^2
      return(fn(w) + pen)
    }

    ui <- diag(k_clean) # w_i >= 0
    ci <- rep(0, k_clean)
    start_w <- rep(1/k_clean, k_clean)

    opt <- try(stats::constrOptim(theta = start_w, f = fn_penalized, grad = NULL,
                                  ui = ui, ci = ci), silent = TRUE)

    if(inherits(opt, "try-error")) {
      final_w <- rep(1/k_clean, k_clean)
    } else {
      final_w <- opt$par
      final_w[final_w < 0] <- 0
      if(sum(final_w) > 0) final_w <- final_w / sum(final_w) else final_w <- rep(1/k_clean, k_clean)
    }
  }

  # Expand back to full vector
  full_w <- rep(0, ncol(preds))
  full_w[valid_cols] <- final_w
  return(full_w)
}




