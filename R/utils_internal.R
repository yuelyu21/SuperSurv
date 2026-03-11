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
                                           newdata = subset(tempValid, select = uniqueScreen[j,], drop = FALSE),
                                           new.times = t.grid,
                                           id = tempId,
                                           obsWeights = tempObsWeights)))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", predAlg,
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      } else {
        libraryRows <- which(library$library$predAlgorithm == predAlg & library$library$rowScreen %in% unlist(screenMap[j]))
        for (row in libraryRows) {
          # --- Add this safety check ---
          expected_length <- nrow(out[, , row]) * ncol(out[, , row])
          actual_length <- length(testAlg$pred)

          if (actual_length != expected_length) {
            stop(sprintf("\n>>> CRASH CAUGHT! <<<\nAlgorithm Index: %s\nExpected matrix size: %s x %s\nBut the algorithm returned length: %s\nCheck if this wrapper was updated to accept 'newdata'!",
                         row, nrow(out[, , row]), ncol(out[, , row]), actual_length))
          }

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
  testScreen <- try(do.call(screen_fn, list), silent = TRUE)

  if (inherits(testScreen, "try-error")) {
    # Print the ACTUAL error so the developer can see it!
    warning("Screening crashed in: ", fun, call. = FALSE)
    warning("Error log: ", testScreen[1], call. = FALSE)
    warning("Defaulting to keeping all variables.", call. = FALSE)
    out <- rep(TRUE, ncol(list$X))
  } else {
    out <- testScreen
  }
  return(out)
}


.predFun <- function(index, lib, time, event, dataX, newdata, whichScreen, t.grid,
                     family, id, obsWeights, verbose, control, libraryNames) {
  if (verbose) {
    message(paste("full", libraryNames[index]))
  }
  pred_fn <- get(lib$predAlgorithm[index])
  testAlg <- try(do.call(pred_fn, list(time = time, event = event,
                                       X = subset(dataX, select = whichScreen[lib$rowScreen[index], ],
                                                  drop = FALSE),
                                       newdata = subset(newdata, select = whichScreen[lib$rowScreen[index],
                                       ], drop = FALSE), id = id,
                                       obsWeights = obsWeights, new.times = t.grid)))
  if (inherits(testAlg, "try-error")) {
    warning(paste("Error in algorithm", lib$predAlgorithm[index],
                  " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
    out <- rep.int(NA, times = nrow(newdata))
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
.checkInputs <- function(time, event, X, newdata, id, obsWeights, verbose) {
  if(any(time < 0)) stop("Only non-negative event/censoring times allowed!")
  if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
  if(any(is.na(time)) | any(is.na(event))) stop("No missing values allowed in time or event.")
  if(any(is.na(X)) | any(is.na(newdata))) stop("No missing values allowed in X or new X.")
  if(length(time) != length(event) | length(time) != nrow(X)) stop("time and event must be n x 1 vectors and X must have n rows.")
  if(!is.data.frame(X) | !is.data.frame(newdata)) stop("X and newdata must be data frames.")
  if(!identical(names(X), names(newdata))) stop("X and newdata must have the same features.")
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
        SL.library[[ii]] <- c(SL.library[[ii]], "")
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











#' Internal function to perform iterative Super Learner optimization
#' @noRd
.surviterativeSL <- function(event.Z, cens.Z, time, event, X, obsWeights, id, control, verbose,
                             event.errorsInLibrary, cens.errorsInLibrary,
                             metalearner = "brier") {

  if (verbose) message("Performing iterative SuperLearner optimization...")
  if (!metalearner %in% c("brier", "entropy", "logloss")) {
    stop("metalearner must be one of: 'brier', 'entropy', 'logloss'")
  }

  trace <- isTRUE(control$traceIter)   # don't depend on verbose
  tol   <- if (!is.null(control$tol)) control$tol else 1e-5

  event.k <- dim(event.Z)[3]
  cens.k  <- dim(cens.Z)[3]
  N <- length(time)
  event.n.time <- length(control$event.t.grid)
  cens.n.time  <- length(control$cens.t.grid)

  epsilon <- max(min(diff(sort(unique(time)))), 1e-5)

  # ---- flatten to long matrices: (N*T) x K ----
  event.Z.long <- matrix(NA_real_, nrow = N * event.n.time, ncol = event.k)
  for (j in seq_len(event.k)) event.Z.long[, j] <- c(event.Z[,,j])

  cens.Z.long <- matrix(NA_real_, nrow = N * cens.n.time, ncol = cens.k)
  for (j in seq_len(cens.k)) cens.Z.long[, j] <- c(cens.Z[,,j])

  # ---- observed values at each subject time (for updating the other side) ----
  event.Z.obs <- matrix(NA_real_, nrow = N, ncol = event.k)
  for (i in seq_len(N)) {
    for (j in seq_len(event.k)) {
      event.Z.obs[i, j] <- stats::approx(
        x = control$event.t.grid, y = event.Z[i, , j],
        xout = time[i], method = "constant", rule = 2, ties = mean
      )$y
    }
  }

  cens.Z.obs <- matrix(NA_real_, nrow = N, ncol = cens.k)
  for (i in seq_len(N)) {
    for (j in seq_len(cens.k)) {
      cens.Z.obs[i, j] <- stats::approx(
        x = c(-1, control$cens.t.grid), y = c(1, cens.Z[i, , j]),
        xout = time[i] - epsilon, method = "constant", rule = 2, ties = mean
      )$y
    }
  }

  # ---- long vectors ----
  obsWeights.event.long <- rep(obsWeights, event.n.time)
  obsWeights.cens.long  <- rep(obsWeights, cens.n.time)
  time.event.long <- rep(time, event.n.time)
  time.cens.long  <- rep(time, cens.n.time)
  event.event.long <- rep(event, event.n.time)
  event.cens.long  <- rep(event, cens.n.time)
  event.t.grid.long <- rep(control$event.t.grid, each = N)
  cens.t.grid.long  <- rep(control$cens.t.grid,  each = N)

  initWeightAlg <- get(control$initWeightAlg)

  obs.cens.vals <- NULL
  obs.event.vals <- NULL
  S.coef <- rep(0, event.k)
  G.coef <- rep(0, cens.k)

  # ---- Initialization ----
  if (control$initWeight == "censoring") {

    initFit <- initWeightAlg(
      time = time, event = 1 - event, X = X, newdata = X,
      new.times = time - epsilon, obsWeights = obsWeights, id = id
    )

    obs.cens.vals <- rep(diag(initFit$pred), length(control$event.t.grid))
    obs.cens.vals <- pmax(obs.cens.vals, 1e-4)

    # logloss needs G(t|X) on the grid which we don't have yet -> init with brier
    init_method <- if (metalearner == "logloss") "brier" else metalearner

    S.coef[!event.errorsInLibrary] <- .survcomputeCoef(
      time = time.event.long, event = event.event.long,
      t.vals = event.t.grid.long, cens.vals = obs.cens.vals,
      preds = event.Z.long[, !event.errorsInLibrary, drop = FALSE],
      obsWeights = obsWeights.event.long,
      method = init_method
    )

    obs.event.vals <- rep(drop(event.Z.obs %*% S.coef), length(control$cens.t.grid))
    obs.event.vals <- pmax(obs.event.vals, 1e-4)

  } else {

    initFit <- initWeightAlg(
      time = time, event = event, X = X, newdata = X,
      new.times = time, obsWeights = obsWeights, id = id
    )

    obs.event.vals <- rep(diag(initFit$pred), length(control$cens.t.grid))
    obs.event.vals <- pmax(obs.event.vals, 1e-4)
  }



  # ---- Iteration ----
  iter <- 1
  G_t_long <- NULL

  while (TRUE) {


    if (iter > control$max.SL.iter) {
      warning("Did not converge in ", control$max.SL.iter, " iterations")
      break
    }

    if (!is.null(obs.cens.vals)) obs.cens.vals.old <- obs.cens.vals
    if (!is.null(obs.event.vals)) obs.event.vals.old <- obs.event.vals

    # ---- Update censoring weights (G.coef) ----
    G.coef <- rep(0, cens.k)
    # IMPORTANT: Do NOT use logloss for censoring weight update.
    # Keep censoring update on Brier (Westling-style).
    G_method <- if (metalearner == "logloss") "brier" else metalearner

    G.coef[!cens.errorsInLibrary] <- .survcomputeCoef(
      time = time.cens.long, event = 1 - event.cens.long,
      t.vals = cens.t.grid.long, cens.vals = obs.event.vals,
      preds = cens.Z.long[, !cens.errorsInLibrary, drop=FALSE],
      obsWeights = obsWeights.cens.long,
      method = G_method
    )

    # ---- Update obs.cens.vals = G(tildeT|X) replicated across event grid ----
    obs.cens.vals <- rep(drop(cens.Z.obs %*% G.coef), length(control$event.t.grid))
    obs.cens.vals <- pmax(obs.cens.vals, 1e-4)

    # ---- If logloss: compute G_t_long = G(t|X) for each (i,t_event) row ----
    if (metalearner == "logloss") {

      if (all(cens.errorsInLibrary)) {
        stop("All censoring learners failed; cannot compute G(t|X) for logloss.")
      }

      # 1) ensemble on censoring grid using LONG matrix (2D) -> no conformability error
      G_cens_long <- drop(
        cens.Z.long[, !cens.errorsInLibrary, drop = FALSE] %*% G.coef[!cens.errorsInLibrary]
      )
      G_cens_long <- pmax(G_cens_long, 1e-4)

      # 2) reshape to N x cens.n.time
      G_cens_grid <- matrix(G_cens_long, nrow = N, ncol = cens.n.time)

      # 3) interpolate each row from censoring grid -> event grid, then flatten long
      G_event_grid <- matrix(NA_real_, nrow = N, ncol = event.n.time)
      for (i in seq_len(N)) {
        G_event_grid[i, ] <- stats::approx(
          x = control$cens.t.grid,
          y = G_cens_grid[i, ],
          xout = control$event.t.grid,
          method = "constant",
          rule = 2,
          ties = mean
        )$y
      }
      G_event_grid <- pmax(G_event_grid, 1e-4)
      G_t_long <- as.vector(G_event_grid)

      if (length(G_t_long) != length(event.t.grid.long)) {
        stop("Internal error: G_t_long length mismatch with event.t.grid.long.")
      }
    }

    # ---- Update event weights (S.coef) ----
    if (metalearner == "logloss") {
      if (is.null(G_t_long)) stop("G_t_long missing for logloss; check censoring library.")

      S.coef <- rep(0, event.k)
      S.coef[!event.errorsInLibrary] <- .survcomputeCoef(
        time = time.event.long, event = event.event.long,
        t.vals = event.t.grid.long,
        cens.vals = obs.cens.vals,      # G(tildeT|X) replicated across t
        G_t = G_t_long,                 # G(t|X) for each (i,t) row
        preds = event.Z.long[, !event.errorsInLibrary, drop = FALSE],
        obsWeights = obsWeights.event.long,
        method = "logloss"
      )

    } else {

      S.coef <- rep(0, event.k)
      S.coef[!event.errorsInLibrary] <- .survcomputeCoef(
        time = time.event.long, event = event.event.long,
        t.vals = event.t.grid.long, cens.vals = obs.cens.vals,
        preds = event.Z.long[, !event.errorsInLibrary, drop = FALSE],
        obsWeights = obsWeights.event.long,
        method = metalearner
      )
    }

    obs.event.vals <- rep(drop(event.Z.obs %*% S.coef), length(control$cens.t.grid))
    obs.event.vals <- pmax(obs.event.vals, 1e-4)



    if (verbose) {
      topK_S <- min(3, length(S.coef))
      topK_G <- min(3, length(G.coef))

      ordS <- order(S.coef, decreasing = TRUE)[seq_len(topK_S)]
      ordG <- order(G.coef, decreasing = TRUE)[seq_len(topK_G)]

      labS <- if (!is.null(names(S.coef))) names(S.coef)[ordS] else as.character(ordS)
      labG <- if (!is.null(names(G.coef))) names(G.coef)[ordG] else as.character(ordG)

      if (isTRUE(verbose)) {
        message(sprintf("[iter=%d] top S: %s | top G: %s",
                  iter,
                  paste(paste0(labS, "=", round(S.coef[ordS], 3)), collapse = ", "),
                  paste(paste0(labG, "=", round(G.coef[ordG], 3)), collapse = ", ")))
      }
      utils::flush.console()
    }

    # ---- Convergence ----
    # ---- after updating obs.cens.vals and obs.event.vals ----
    cens.delta  <- Inf
    event.delta <- Inf

    if (!is.null(obs.cens.vals.old) && !is.null(obs.event.vals.old)) {
      cens.delta  <- max(abs(obs.cens.vals  - obs.cens.vals.old))
      event.delta <- max(abs(obs.event.vals - obs.event.vals.old))
    }

    if (isTRUE(verbose)) {
      message(sprintf("[iter=%d] cens.delta=%.3e  event.delta=%.3e  sum=%.3e",
                iter, cens.delta, event.delta, cens.delta + event.delta))
      }

    tol <- if (!is.null(control$tol)) control$tol else 1e-5
    if (cens.delta + event.delta < tol) {
      if (verbose) message("Converged in ", iter, " iterations.")
      break
    }



    iter <- iter + 1
  }

  # ---- CV risks ----
  calc_risk <- function(preds, time, event, t.grid, cens_T, weights, method, G_t = NULL) {

    if (method == "brier") {
      Y_ipcw <- 1 - (as.numeric(time <= t.grid) * event / cens_T)
      return(apply(preds, 2, function(col) mean(weights * (Y_ipcw - col)^2, na.rm = TRUE)))
    }

    if (method == "entropy") {
      Y_ipcw <- 1 - (as.numeric(time <= t.grid) * event / cens_T)
      Y_clamp <- pmax(pmin(Y_ipcw, 1), 0)
      return(apply(preds, 2, function(col) {
        p <- pmax(pmin(col, 1 - 1e-15), 1e-15)
        loss <- -(weights * (Y_clamp * log(p) + (1 - Y_clamp) * log(1 - p)))
        mean(loss, na.rm = TRUE)
      }))
    }

    if (method == "logloss") {
      if (is.null(G_t)) stop("calc_risk(method='logloss') requires G_t.")
      G_T <- pmax(cens_T, 1e-4)
      G_t <- pmax(G_t, 1e-4)
      fail <- as.numeric(time <= t.grid) * event
      surv <- as.numeric(time >  t.grid)

      w_cap <- 100
      w_fail <- pmin(fail / G_T, w_cap)
      w_surv <- pmin(surv / G_t, w_cap)

      return(apply(preds, 2, function(col) {
        S <- pmin(pmax(col, 1e-10), 1 - 1e-10)
        loss <- -weights * (w_fail * log(1 - S) + w_surv * log(S))
        mean(loss, na.rm = TRUE)
      }))
    }

    stop("Unknown method in calc_risk")
  }

  # for risk calc, rebuild final G_t_long the same way if logloss
  if (metalearner == "logloss") {
    if (is.null(G_t_long)) stop("G_t_long missing at risk computation stage.")
    G_t_long_final <- G_t_long
  } else {
    G_t_long_final <- NULL
  }

  event.cvRisks <- calc_risk(
    preds = event.Z.long,
    time = time.event.long,
    event = event.event.long,
    t.grid = event.t.grid.long,
    cens_T = obs.cens.vals,
    weights = obsWeights.event.long,
    method = metalearner,
    G_t = G_t_long_final
  )

  # keep censoring risk as brier (same as Westling-style)
  cens.cvRisks <- calc_risk(
    preds = cens.Z.long,
    time = time.cens.long,
    event = event.cens.long,
    t.grid = cens.t.grid.long,
    cens_T = obs.event.vals,
    weights = obsWeights.cens.long,
    method = "brier"
  )

  list(
    event.coef = S.coef,
    cens.coef  = G.coef,
    event.cvRisks = event.cvRisks,
    cens.cvRisks  = cens.cvRisks
  )
}











#' Internal function to compute ensemble weights
#' @importFrom nnls nnls
#' @noRd
.survcomputeCoef <- function(time, event, t.vals, cens.vals, preds, obsWeights, method = "brier", G_t = NULL) {

  if (!method %in% c("brier","entropy","logloss")) {
    stop("Unknown method: ", method, ". Use 'brier', 'entropy', or 'logloss'.")
  }

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
  if (method == "brier") {
    # IPCW "Observed" Outcome (1 if alive, 0 if dead, weighted)
    out <- 1 - as.numeric(time <= t.vals) * event / cens.vals

    fit.nnls <- nnls::nnls(sqrt(obsWeights) * preds_clean, sqrt(obsWeights) * out)
    coef <- coef(fit.nnls)
    if(sum(coef) == 0) coef <- rep(1, length(coef))

    final_w <- coef / sum(coef)
  }


  # ----------------------------------------------------------------------------
  # METHOD 2: entropy (cross-entropy on IPCW pseudo-outcomes; simplex via softmax)
  # ----------------------------------------------------------------------------
  else if (method == "entropy") {

    # Raw IPCW pseudo-outcome
    Y_ipcw <- 1 - (as.numeric(time <= t.vals) * event / cens.vals)

    # Clamp target into [0,1] so cross-entropy is well-defined
    Y_clamp <- pmax(pmin(Y_ipcw, 1), 0)

    # Softmax maps unconstrained theta -> simplex weights
    softmax <- function(theta) {
      z <- theta - max(theta)
      w <- exp(z)
      w / sum(w)
    }

    # Loss as a function of weights w (same as your original idea)
    loss_w <- function(w) {
      P_ens <- drop(preds_clean %*% w)
      P_ens <- pmax(pmin(P_ens, 1 - 1e-10), 1e-10)

      loss <- -(obsWeights * (Y_clamp * log(P_ens) + (1 - Y_clamp) * log(1 - P_ens)))
      mean(loss)
    }

    # Optimize over theta (unconstrained), then transform via softmax
    loss_theta <- function(theta) loss_w(softmax(theta))

    opt <- try(stats::optim(par = rep(0, k_clean), fn = loss_theta, method = "BFGS"),
               silent = TRUE)

    if (inherits(opt, "try-error")) {
      final_w <- rep(1 / k_clean, k_clean)
    } else {
      final_w <- softmax(opt$par)
    }
  }

  # ----------------------------------------------------------------------------
  # METHOD 3: IPCW log-loss (simplex via softmax + optim)
  # ----------------------------------------------------------------------------
  else if (method == "logloss") {

    if (is.null(G_t)) {
      stop("method='logloss' requires G_t = G(t|X) aligned with the (i,t) long rows.")
    }

    # Floors for numerical stability
    G_T <- pmax(cens.vals, 1e-4)   # should be G(tildeT|X) replicated across t
    G_t <- pmax(G_t, 1e-4)         # should be G(t|X) for each (i,t) row

    fail <- as.numeric(time <= t.vals) * event      # Δ I(T<=t)
    surv <- as.numeric(time >  t.vals)              # I(T>t)

    # Softmax maps unconstrained theta -> simplex weights
    softmax <- function(theta) {
      z <- theta - max(theta)
      w <- exp(z)
      w / sum(w)
    }

    # Loss as a function of weights w
    loss_w <- function(w) {
      S_ens <- drop(preds_clean %*% w)
      S_ens <- pmin(pmax(S_ens, 1e-10), 1 - 1e-10)

      # Optional but recommended: cap IPCW multipliers
      w_cap <- 100
      w_fail <- pmin(fail / G_T, w_cap)
      w_surv <- pmin(surv / G_t, w_cap)

      loss <- -obsWeights * (w_fail * log(1 - S_ens) + w_surv * log(S_ens))
      mean(loss)
    }

    # Optimize over theta (unconstrained), then transform via softmax
    loss_theta <- function(theta) loss_w(softmax(theta))

    opt <- try(stats::optim(par = rep(0, k_clean), fn = loss_theta, method = "BFGS"),
               silent = TRUE)

    if (inherits(opt, "try-error")) {
      final_w <- rep(1 / k_clean, k_clean)
    } else {
      final_w <- softmax(opt$par)
    }
  }

  # Expand back to full vector
  full_w <- rep(0, ncol(preds))
  full_w[valid_cols] <- final_w
  return(full_w)
}




