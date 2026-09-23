########## function previously inside the main function
### Get cross-validated survivals for S and T
#' @noRd
.crossValFUN <- function(valid, validRows, time, event, dataX, id, obsWeights, t.grid,
                         library, kScreen, k, p, verbose, functions) {

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
    screen_fn <- functions[[library$screenAlgorithm[s]]]
    testScreen <- try(do.call(screen_fn, list(time = tempTime,
                                              event = tempEvent,
                                              X = tempLearn, id = tempId,
                                              obsWeights = tempObsWeights)))
    if (inherits(testScreen, "try-error")) {
      warning(
        "Screener '", library$screenAlgorithm[s], "' failed in cross-validation fold ",
        foldNum, ": ", .try_error_message(testScreen),
        ". All variables will be retained for this fold.",
        call. = FALSE
      )
      tempWhichScreen[s, ] <- TRUE
    }
    else {
      tempWhichScreen[s, ] <- .validate_screener_output(
        testScreen,
        screener = library$screenAlgorithm[s],
        n_features = p,
        context = paste0("cross-validation fold ", foldNum)
      )
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
    pred_fn <- functions[[predAlg]]
    for(j in seq(nrow(uniqueScreen))) {
      testAlg <- try(do.call(pred_fn, list(time = tempTime, event = tempEvent,
                                           X = subset(tempLearn, select = uniqueScreen[j,], drop = FALSE),
                                           newdata = subset(tempValid, select = uniqueScreen[j,], drop = FALSE),
                                           new.times = t.grid,
                                           id = tempId,
                                           obsWeights = tempObsWeights)))
      if (inherits(testAlg, "try-error")) {
        warning(
          "Learner '", predAlg, "' failed in cross-validation fold ", foldNum,
          ": ", .try_error_message(testAlg),
          ". It will receive zero ensemble weight.",
          call. = FALSE
        )
      } else {
        learner_pred <- .validate_learner_output(
          testAlg,
          learner = predAlg,
          n_observations = nrow(tempValid),
          times = t.grid,
          context = paste0("cross-validation fold ", foldNum)
        )
        libraryRows <- which(library$library$predAlgorithm == predAlg & library$library$rowScreen %in% unlist(screenMap[j]))
        for (row in libraryRows) {
          out[, , row] <- learner_pred
        }
      }
    }
  }


  invisible(list(out = out))
}




# Generate Full Screening logic
#' @noRd
.screenFun <- function(fun, list, functions) {
  screen_fn <- functions[[fun]]
  # Robust screen call
  testScreen <- try(do.call(screen_fn, list), silent = TRUE)

  if (inherits(testScreen, "try-error")) {
    warning(
      "Screener '", fun, "' failed on the full training data: ",
      .try_error_message(testScreen), ". All variables will be retained.",
      call. = FALSE
    )
    out <- rep(TRUE, ncol(list$X))
  } else {
    out <- .validate_screener_output(
      testScreen,
      screener = fun,
      n_features = ncol(list$X),
      context = "the full training data"
    )
  }
  return(out)
}


.predFun <- function(index, lib, time, event, dataX, newdata, whichScreen, t.grid,
                     family, id, obsWeights, verbose, control, libraryNames, functions) {
  if (verbose) {
    message(paste("full", libraryNames[index]))
  }
  pred_fn <- functions[[lib$predAlgorithm[index]]]
  testAlg <- try(do.call(pred_fn, list(time = time, event = event,
                                       X = subset(dataX, select = whichScreen[lib$rowScreen[index], ],
                                                  drop = FALSE),
                                       newdata = subset(newdata, select = whichScreen[lib$rowScreen[index],
                                       ], drop = FALSE), id = id,
                                       obsWeights = obsWeights, new.times = t.grid)))
  if (inherits(testAlg, "try-error")) {
    warning(
      "Learner '", lib$predAlgorithm[index], "' failed on the full training data: ",
      .try_error_message(testAlg), ". It will receive zero ensemble weight.",
      call. = FALSE
    )
    out <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(t.grid))
    model_out <- NULL
  }
  else {
    out <- .validate_learner_output(
      testAlg,
      learner = lib$predAlgorithm[index],
      n_observations = nrow(newdata),
      times = t.grid,
      require_fit = isTRUE(control$saveFitLibrary),
      context = "the full training data"
    )
    if(control$saveFitLibrary) model_out <- testAlg$fit
    else model_out <- NULL
  }

  invisible(list(out = out, model_out = model_out))
}






########## function previously inside the main function
#' @noRd
.checkInputs <- function(time, event, X, newdata, id, obsWeights, verbose) {
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(time) || length(time) == 0L || any(!is.finite(time))) {
    stop("`time` must be a non-empty numeric vector containing only finite values.", call. = FALSE)
  }
  if (any(time < 0)) {
    stop("`time` must contain only non-negative values.", call. = FALSE)
  }
  if (!is.numeric(event) || length(event) != length(time) ||
      any(!is.finite(event)) || any(!event %in% c(0, 1))) {
    stop("`event` must be a numeric 0/1 vector with the same length as `time`.", call. = FALSE)
  }
  if (!any(event == 1)) {
    stop("`event` must contain at least one observed event (1).", call. = FALSE)
  }
  if (!any(event == 0)) {
    stop("`event` must contain at least one censored observation (0).", call. = FALSE)
  }
  if (nrow(X) != length(time)) {
    stop("`X` must have one row for each value in `time` and `event`.", call. = FALSE)
  }
  if (nrow(X) == 0L || ncol(X) == 0L) {
    stop("`X` must contain at least one row and one feature column.", call. = FALSE)
  }
  .validate_data_frame_columns(X, "X")
  .validate_data_frame_columns(newdata, "newdata")
  if (nrow(newdata) == 0L) {
    stop("`newdata` must contain at least one row.", call. = FALSE)
  }
  if (!identical(names(X), names(newdata))) {
    stop("`newdata` must contain exactly the same feature columns as `X`, in the same order.",
         call. = FALSE)
  }
  if (!is.numeric(obsWeights) || length(obsWeights) != length(time) ||
      any(!is.finite(obsWeights)) || any(obsWeights < 0) || sum(obsWeights) <= 0) {
    stop("`obsWeights` must be a finite, non-negative numeric vector with one value per observation and a positive sum.",
         call. = FALSE)
  }
  .validate_scalar_logical(verbose, "verbose")
  if (!is.null(id) && length(id) != length(time)) {
    stop("`id` must be NULL or have the same length as `time`.", call. = FALSE)
  }
}


#' @noRd
.validate_scalar_logical <- function(value, arg_name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", arg_name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}


#' @noRd
.validate_time_grid <- function(times, arg_name = "times") {
  if (!is.numeric(times) || length(times) == 0L || any(!is.finite(times))) {
    stop("`", arg_name, "` must be a non-empty numeric vector containing only finite values.",
         call. = FALSE)
  }
  if (any(times < 0)) {
    stop("`", arg_name, "` must contain only non-negative values.", call. = FALSE)
  }
  if (is.unsorted(times, strictly = TRUE)) {
    stop("`", arg_name, "` must be strictly increasing with no duplicate values.",
         call. = FALSE)
  }
  as.numeric(times)
}


#' @noRd
.validate_data_frame_columns <- function(data, arg_name) {
  if (is.null(names(data)) || anyNA(names(data)) || any(!nzchar(names(data))) ||
      anyDuplicated(names(data))) {
    stop("`", arg_name, "` must have unique, non-empty column names.", call. = FALSE)
  }
  if (anyNA(data)) {
    stop("`", arg_name, "` must not contain missing values.", call. = FALSE)
  }
  numeric_columns <- vapply(data, is.numeric, logical(1L))
  if (any(numeric_columns) &&
      any(!vapply(data[numeric_columns], function(x) all(is.finite(x)), logical(1L)))) {
    stop("Numeric columns in `", arg_name, "` must contain only finite values.",
         call. = FALSE)
  }
  invisible(TRUE)
}


#' @noRd
.validate_SuperSurv_object <- function(object, arg_name = "object") {
  if (!inherits(object, "SuperSurv") || !is.list(object)) {
    stop("`", arg_name, "` must be a fitted 'SuperSurv' object returned by `SuperSurv()`.",
         call. = FALSE)
  }
  required <- c(
    "event.coef", "cens.coef", "event.libraryNames", "cens.libraryNames",
    "control", "varNames"
  )
  missing_components <- required[!vapply(required, function(name) {
    !is.null(object[[name]])
  }, logical(1L))]
  if (length(missing_components)) {
    stop(
      "`", arg_name, "` is not a complete fitted 'SuperSurv' object; missing component(s): ",
      paste(missing_components, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' @noRd
.validate_prediction_newdata <- function(newdata, object) {
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  if (nrow(newdata) == 0L) {
    stop("`newdata` must contain at least one row.", call. = FALSE)
  }
  .validate_data_frame_columns(newdata, "newdata")
  required <- object$varNames
  missing_features <- setdiff(required, names(newdata))
  unexpected_features <- setdiff(names(newdata), required)
  if (length(missing_features)) {
    stop("`newdata` is missing training feature(s): ",
         paste(missing_features, collapse = ", "), ".", call. = FALSE)
  }
  if (length(unexpected_features)) {
    stop("`newdata` contains unexpected feature(s): ",
         paste(unexpected_features, collapse = ", "),
         ". Supply exactly the variables used to fit the model.", call. = FALSE)
  }
  newdata[, required, drop = FALSE]
}


#' @noRd
.validate_learner_output <- function(result, learner, n_observations, times,
                                     require_fit = FALSE, context = "fitting") {
  prefix <- paste0("Learner '", learner, "' returned an invalid result during ", context, ": ")
  if (!is.list(result)) {
    stop(prefix, "expected a list with components `pred` and `fit`.", call. = FALSE)
  }
  if (is.null(result$pred)) {
    stop(prefix, "missing required component `pred`.", call. = FALSE)
  }
  if (!is.matrix(result$pred) || !is.numeric(result$pred)) {
    stop(prefix, "`pred` must be a numeric matrix.", call. = FALSE)
  }
  expected <- c(as.integer(n_observations), as.integer(length(times)))
  if (!identical(dim(result$pred), expected)) {
    stop(
      prefix, "`pred` has dimensions ", paste(dim(result$pred), collapse = " x "),
      "; expected ", paste(expected, collapse = " x "),
      " (observations x prediction times).",
      call. = FALSE
    )
  }
  if (any(!is.finite(result$pred))) {
    stop(prefix, "`pred` must contain only finite values.", call. = FALSE)
  }
  if (any(result$pred < 0 | result$pred > 1)) {
    observed <- range(result$pred)
    stop(
      prefix, "`pred` must contain survival probabilities in [0, 1]; observed range was [",
      format(observed[1L]), ", ", format(observed[2L]), "].",
      call. = FALSE
    )
  }
  if (isTRUE(require_fit) && is.null(result$fit)) {
    stop(prefix, "missing required component `fit` while `saveFitLibrary = TRUE`.",
         call. = FALSE)
  }
  result$pred
}


#' @noRd
.validate_screener_output <- function(result, screener, n_features,
                                      context = "fitting") {
  prefix <- paste0("Screener '", screener, "' returned an invalid result during ", context, ": ")
  if (!is.logical(result) || length(result) != n_features || anyNA(result)) {
    stop(
      prefix, "expected a complete logical vector of length ", n_features,
      ", one value per feature.",
      call. = FALSE
    )
  }
  if (!any(result)) {
    stop(prefix, "at least one feature must be selected.", call. = FALSE)
  }
  result
}


#' @noRd
.try_error_message <- function(error) {
  message <- as.character(error)[1L]
  sub("^Error[^:]*:\\s*", "", message)
}


#' @noRd
.find_library_function <- function(name, envir) {
  fun <- get0(name, envir = envir, mode = "function", inherits = TRUE)
  if (is.null(fun)) {
    fun <- get0(name, envir = environment(.find_library_function),
                mode = "function", inherits = TRUE)
  }
  fun
}


#' @noRd
.resolve_library_functions <- function(library, envir) {
  names <- unique(c(library$library$predAlgorithm, library$screenAlgorithm))
  stats::setNames(lapply(names, .find_library_function, envir = envir), names)
}


#' @noRd
.validate_library_argument <- function(survSL.library, arg_name, envir = parent.frame()) {
  format_spec <- paste0(
    "`", arg_name, "` must be either a character vector of learner names ",
    "(for example, `c(\"surv.coxph\", \"surv.km\")`) or a list where each element ",
    "is a character vector beginning with a learner name followed optionally by screener names ",
    "(for example, `list(c(\"surv.coxph\", \"screen.all\"))`)."
  )

  validate_names <- function(x, context) {
    context_label <- if (nzchar(context)) paste0(context, " ") else ""
    if (!is.character(x)) {
      stop("`", arg_name, "` ", context_label, "must be a character vector. ", format_spec, call. = FALSE)
    }
    if (length(x) < 1L) {
      stop("`", arg_name, "` ", context_label, "must contain at least one function name. ", format_spec, call. = FALSE)
    }
    if (any(is.na(x)) || any(!nzchar(trimws(x)))) {
      stop("`", arg_name, "` ", context_label, "contains missing or empty function names. ", format_spec, call. = FALSE)
    }

    missing_funs <- x[vapply(x, function(name) {
      is.null(.find_library_function(name, envir))
    }, logical(1L))]
    if (length(missing_funs) > 0L) {
      stop(
        "`", arg_name, "` ", context_label, "contains unknown function name(s): ",
        paste(unique(missing_funs), collapse = ", "),
        ". Define these functions before calling `SuperSurv()`. ",
        format_spec,
        call. = FALSE
      )
    }
  }

  if (inherits(survSL.library, "SuperSurv_grid")) {
    validate_names(as.character(survSL.library), "generated by `create_grid()`")
    return(invisible(TRUE))
  }

  if (is.character(survSL.library)) {
    validate_names(survSL.library, "")
    return(invisible(TRUE))
  }

  if (!is.list(survSL.library)) {
    stop(format_spec, call. = FALSE)
  }

  if (length(survSL.library) < 1L) {
    stop("`", arg_name, "` must not be empty. ", format_spec, call. = FALSE)
  }

  for (ii in seq_along(survSL.library)) {
    validate_names(survSL.library[[ii]], paste0("element ", ii))
  }

  invisible(TRUE)
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
        survSL.library[[ii]] <- c(survSL.library[[ii]], "screen.all")
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
.surviterativeSL_legacy <- function(event.Z, cens.Z, time, event, X, obsWeights, id, control, verbose,
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
.survcomputeCoef_legacy <- function(time, event, t.vals, cens.vals, preds, obsWeights, method = "brier", G_t = NULL) {

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


#' Evaluate a right-continuous step curve at one or more times
#' @noRd
.step_curve_at <- function(grid, values, xout, left.limit = FALSE, initial = 1) {
  if (length(grid) != length(values)) {
    stop("Internal error: prediction grid and curve have different lengths.")
  }
  if (length(grid) == 0L) return(rep(initial, length(xout)))

  index <- vapply(xout, function(x) {
    locations <- which(if (left.limit) grid < x else grid <= x)
    if (length(locations)) max(locations) else 0L
  }, integer(1L))

  out <- rep(initial, length(index))
  has_value <- index > 0L
  out[has_value] <- values[index[has_value]]
  out
}


#' Evaluate subject-specific step curves at subject-specific times
#' @noRd
.row_step_curve_at <- function(pred, grid, xout, left.limit = FALSE,
                               initial = 1) {
  pred <- as.matrix(pred)
  if (nrow(pred) != length(xout) || ncol(pred) != length(grid)) {
    stop("Internal error: prediction matrix is not aligned with subjects and times.")
  }
  vapply(seq_along(xout), function(i) {
    .step_curve_at(grid, pred[i, ], xout[i], left.limit, initial)
  }, numeric(1L))
}


#' Stabilize inverse-probability weights and retain diagnostics
#' @noRd
.stabilize_ipcw <- function(probability, floor, cap) {
  if (any(!is.finite(probability))) {
    stop("IPCW denominator probabilities must be finite.")
  }
  floored_probability <- pmax(probability, floor)
  raw_weight <- ifelse(probability > 0, 1 / probability, Inf)
  floored_weight <- 1 / floored_probability
  capped <- floored_weight > cap
  list(
    probability = probability,
    raw.weight = raw_weight,
    weight = pmin(floored_weight, cap),
    floored = probability < floor,
    capped = capped,
    effective.denominator.threshold = max(
      floor, if (is.finite(cap)) 1 / cap else 0
    ),
    floor.fraction = mean(probability < floor),
    cap.fraction = mean(capped)
  )
}


#' Summarize the IPCW contribution used by one loss component
#' @noRd
.summarize_ipcw <- function(stabilized, active, loss.contribution = NULL) {
  active <- as.logical(active)
  if (length(active) != length(stabilized$weight)) {
    stop("Internal error: IPCW diagnostic mask is not aligned with weights.")
  }
  if (!any(active)) {
    return(list(
      rows = 0L, denominator.minimum = NA_real_, floor.count = 0L,
      floor.fraction = 0, cap.count = 0L, cap.fraction = 0,
      raw.weight.quantiles = rep(NA_real_, 6L),
      stabilized.weight.quantiles = rep(NA_real_, 6L),
      effective.sample.size = NA_real_, capped.loss.fraction = 0
    ))
  }

  probabilities <- stabilized$probability[active]
  raw_weights <- stabilized$raw.weight[active]
  weights <- stabilized$weight[active]
  probs <- c(0, 0.5, 0.9, 0.95, 0.99, 1)
  quantile_names <- c("min", "median", "p90", "p95", "p99", "max")
  raw_quantiles <- stats::quantile(raw_weights, probs, names = FALSE, type = 1)
  stabilized_quantiles <- stats::quantile(
    weights, probs, names = FALSE, type = 1
  )
  names(raw_quantiles) <- quantile_names
  names(stabilized_quantiles) <- quantile_names
  weight_sum <- sum(weights)
  effective_sample_size <- if (weight_sum > 0) {
    weight_sum^2 / sum(weights^2)
  } else {
    NA_real_
  }

  capped_loss_fraction <- 0
  if (!is.null(loss.contribution)) {
    contribution <- loss.contribution[active]
    total <- sum(contribution)
    if (is.finite(total) && total > 0) {
      capped_loss_fraction <- sum(
        contribution[stabilized$capped[active]]
      ) / total
    }
  }

  list(
    rows = sum(active),
    denominator.minimum = min(probabilities),
    floor.count = sum(stabilized$floored[active]),
    floor.fraction = mean(stabilized$floored[active]),
    cap.count = sum(stabilized$capped[active]),
    cap.fraction = mean(stabilized$capped[active]),
    raw.weight.quantiles = raw_quantiles,
    stabilized.weight.quantiles = stabilized_quantiles,
    effective.sample.size = effective_sample_size,
    capped.loss.fraction = capped_loss_fraction
  )
}


#' Project a numeric vector onto the probability simplex
#' @noRd
.project_simplex <- function(x) {
  ordered <- sort(x, decreasing = TRUE)
  cumulative <- cumsum(ordered) - 1
  rho <- max(which(ordered - cumulative / seq_along(ordered) > 0))
  theta <- cumulative[rho] / rho
  pmax(x - theta, 0)
}


#' Minimize weighted squared error directly over the probability simplex
#' @noRd
.simplex_least_squares <- function(preds, outcome, weights, tol, maxit) {
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("Observation weights must have a positive finite sum.")
  }

  gram <- crossprod(preds, weights * preds) / weight_sum
  linear <- drop(crossprod(preds, weights * outcome) / weight_sum)
  largest_eigenvalue <- max(eigen(gram, symmetric = TRUE, only.values = TRUE)$values)
  k <- ncol(preds)
  current <- rep(1 / k, k)

  if (!is.finite(largest_eigenvalue) || largest_eigenvalue <= .Machine$double.eps) {
    return(list(
      weights = current, converged = TRUE, iterations = 0L,
      objective = mean(weights * (outcome - drop(preds %*% current))^2)
    ))
  }

  step <- 1 / (2 * largest_eigenvalue)

  if (k == 2L) {
    first_weight <- 0.5
    curvature <- gram[1L, 1L] - 2 * gram[1L, 2L] + gram[2L, 2L]
    linear_term <- linear[1L] - linear[2L] - gram[1L, 2L] + gram[2L, 2L]
    converged <- FALSE
    iterations <- maxit
    for (iteration in seq_len(maxit)) {
      updated <- first_weight - step * (
        curvature * first_weight - linear_term
      )
      updated <- min(max(updated, 0), 1)
      if (abs(updated - first_weight) <= tol) {
        first_weight <- updated
        converged <- TRUE
        iterations <- iteration
        break
      }
      first_weight <- updated
    }
    current <- c(first_weight, 1 - first_weight)
    return(list(
      weights = current,
      converged = converged,
      iterations = iterations,
      objective = mean(weights * (outcome - drop(preds %*% current))^2)
    ))
  }

  converged <- FALSE
  iterations <- maxit
  for (iteration in seq_len(maxit)) {
    gradient <- 2 * drop(gram %*% current - linear)
    updated <- .project_simplex(current - step * gradient)
    if (max(abs(updated - current)) <= tol) {
      current <- updated
      converged <- TRUE
      iterations <- iteration
      break
    }
    current <- updated
  }

  list(
    weights = current,
    converged = converged,
    iterations = iterations,
    objective = mean(weights * (outcome - drop(preds %*% current))^2)
  )
}


#' Optimize a loss over simplex weights using an identifiable softmax map
#' @noRd
.simplex_softmax_optim <- function(k, loss, tol, maxit) {
  if (k == 1L) {
    return(list(weights = 1, converged = TRUE, iterations = 0L,
                objective = loss(1)))
  }
  softmax <- function(theta) {
    logits <- c(theta, 0)
    logits <- logits - max(logits)
    values <- exp(logits)
    values / sum(values)
  }

  opt <- try(
    stats::optim(
      par = rep(0, k - 1L), fn = function(theta) loss(softmax(theta)),
      method = "BFGS", control = list(reltol = tol, maxit = maxit)
    ),
    silent = TRUE
  )
  if (inherits(opt, "try-error") || !is.finite(opt$value)) {
    return(list(weights = rep(1 / k, k), converged = FALSE,
                iterations = NA_integer_, objective = NA_real_))
  }

  list(
    weights = softmax(opt$par),
    converged = identical(opt$convergence, 0L),
    iterations = unname(opt$counts[["function"]]),
    objective = opt$value
  )
}


#' Internal function to compute ensemble weights
#' @noRd
.survcomputeCoef <- function(time, event, t.vals, cens.vals, preds, obsWeights,
                             method = "brier", G_t = NULL,
                             event.inclusive = TRUE, ipcw.floor = 1e-4,
                             ipcw.cap = 100, prob.eps = 1e-10,
                             optimizer.tol = 1e-8,
                             optimizer.maxit = 10000L) {
  if (!method %in% c("brier", "entropy", "logloss")) {
    stop("Unknown method: ", method, ". Use 'brier', 'entropy', or 'logloss'.")
  }

  preds <- as.matrix(preds)
  valid_cols <- !apply(preds, 2, anyNA)
  if (!any(valid_cols)) stop("All candidate prediction columns contain missing values.")
  preds_clean <- preds[, valid_cols, drop = FALSE]
  if (any(!is.finite(preds_clean))) stop("Candidate predictions must be finite.")

  indicator <- if (event.inclusive) time <= t.vals else time < t.vals
  stabilized_T <- .stabilize_ipcw(cens.vals, ipcw.floor, ipcw.cap)
  pseudo_outcome <- 1 - as.numeric(indicator) * event * stabilized_T$weight
  diagnostics <- list(
    method = method, converged = TRUE, iterations = 0L,
    objective = NA_real_, floor.fraction = stabilized_T$floor.fraction,
    cap.fraction = stabilized_T$cap.fraction,
    effective.denominator.threshold =
      stabilized_T$effective.denominator.threshold
  )
  stabilized_t <- NULL

  k_clean <- ncol(preds_clean)
  if (k_clean == 1L && method != "logloss") {
    final_w <- 1
  } else if (method == "brier") {
    fit <- .simplex_least_squares(
      preds_clean, pseudo_outcome, obsWeights, optimizer.tol, optimizer.maxit
    )
    final_w <- fit$weights
    diagnostics[names(fit)[-1L]] <- fit[-1L]
  } else if (method == "entropy") {
    target <- pmax(pmin(pseudo_outcome, 1), 0)
    loss <- function(w) {
      probability <- pmin(pmax(drop(preds_clean %*% w), prob.eps), 1 - prob.eps)
      mean(-obsWeights * (target * log(probability) +
                            (1 - target) * log(1 - probability)))
    }
    fit <- .simplex_softmax_optim(k_clean, loss, optimizer.tol, optimizer.maxit)
    final_w <- fit$weights
    diagnostics[names(fit)[-1L]] <- fit[-1L]
  } else {
    if (is.null(G_t)) {
      stop("method='logloss' requires G_t = G(t|X) aligned with the (i,t) rows.")
    }
    stabilized_t <- .stabilize_ipcw(G_t, ipcw.floor, ipcw.cap)
    fail_weight <- as.numeric(time <= t.vals) * event * stabilized_T$weight
    survive_weight <- as.numeric(time > t.vals) * stabilized_t$weight
    loss <- function(w) {
      survival_probability <- pmin(
        pmax(drop(preds_clean %*% w), prob.eps), 1 - prob.eps
      )
      mean(-obsWeights * (
        fail_weight * log(1 - survival_probability) +
          survive_weight * log(survival_probability)
      ))
    }
    fit <- .simplex_softmax_optim(k_clean, loss, optimizer.tol, optimizer.maxit)
    final_w <- fit$weights
    diagnostics[names(fit)[-1L]] <- fit[-1L]
    diagnostics$floor.fraction <- max(
      diagnostics$floor.fraction, stabilized_t$floor.fraction
    )
    diagnostics$cap.fraction <- max(
      diagnostics$cap.fraction, stabilized_t$cap.fraction
    )
  }

  full_w <- rep(0, ncol(preds))
  full_w[valid_cols] <- final_w
  if (method == "logloss") {
    survival_probability <- pmin(
      pmax(drop(preds_clean %*% final_w), prob.eps), 1 - prob.eps
    )
    failure_active <- time <= t.vals & event == 1
    survivor_active <- time > t.vals
    failure_contribution <- -obsWeights * fail_weight *
      log(1 - survival_probability)
    survivor_contribution <- -obsWeights * survive_weight *
      log(survival_probability)
    diagnostics$failure.weights <- .summarize_ipcw(
      stabilized_T, failure_active, failure_contribution
    )
    diagnostics$survivor.weights <- .summarize_ipcw(
      stabilized_t, survivor_active, survivor_contribution
    )
    combined_contribution <- failure_contribution + survivor_contribution
    combined_active <- failure_active | survivor_active
    combined_probability <- ifelse(
      failure_active, stabilized_T$probability, stabilized_t$probability
    )
    combined_raw_weight <- ifelse(
      failure_active, stabilized_T$raw.weight, stabilized_t$raw.weight
    )
    combined_weight <- ifelse(
      failure_active, stabilized_T$weight, stabilized_t$weight
    )
    combined_floored <- ifelse(
      failure_active, stabilized_T$floored, stabilized_t$floored
    )
    combined_capped <- ifelse(
      failure_active, stabilized_T$capped, stabilized_t$capped
    )
    diagnostics$combined.weights <- .summarize_ipcw(
      list(
        probability = combined_probability,
        raw.weight = combined_raw_weight,
        weight = combined_weight,
        floored = combined_floored,
        capped = combined_capped
      ),
      combined_active,
      combined_contribution
    )
    diagnostics$floor.fraction <- diagnostics$combined.weights$floor.fraction
    diagnostics$cap.fraction <- diagnostics$combined.weights$cap.fraction
  } else {
    diagnostics$weighted.rows <- .summarize_ipcw(
      stabilized_T, as.numeric(indicator) * event > 0
    )
    diagnostics$floor.fraction <- diagnostics$weighted.rows$floor.fraction
    diagnostics$cap.fraction <- diagnostics$weighted.rows$cap.fraction
  }
  attr(full_w, "diagnostics") <- diagnostics
  full_w
}


#' Internal function to perform iterative Super Learner optimization
#' @noRd
.surviterativeSL <- function(event.Z, cens.Z, time, event, X, obsWeights, id,
                             control, verbose, event.errorsInLibrary,
                             cens.errorsInLibrary, metalearner = "brier") {
  if (verbose) message("Performing iterative SuperLearner optimization...")
  if (!metalearner %in% c("brier", "entropy", "logloss")) {
    stop("metalearner must be one of: 'brier', 'entropy', 'logloss'")
  }
  if (all(event.errorsInLibrary)) stop("All event learners failed during cross-validation.")
  if (all(cens.errorsInLibrary)) stop("All censoring learners failed during cross-validation.")

  event.k <- dim(event.Z)[3]
  cens.k <- dim(cens.Z)[3]
  N <- length(time)
  event.n.time <- length(control$event.t.grid)
  cens.n.time <- length(control$cens.t.grid)

  event.Z.long <- matrix(NA_real_, N * event.n.time, event.k)
  for (j in seq_len(event.k)) event.Z.long[, j] <- c(event.Z[, , j])
  cens.Z.long <- matrix(NA_real_, N * cens.n.time, cens.k)
  for (j in seq_len(cens.k)) cens.Z.long[, j] <- c(cens.Z[, , j])

  event.Z.obs <- matrix(NA_real_, N, event.k)
  for (i in seq_len(N)) {
    for (j in seq_len(event.k)) {
      event.Z.obs[i, j] <- .step_curve_at(
        control$event.t.grid, event.Z[i, , j], time[i]
      )
    }
  }
  cens.Z.obs <- matrix(NA_real_, N, cens.k)
  for (i in seq_len(N)) {
    for (j in seq_len(cens.k)) {
      cens.Z.obs[i, j] <- .step_curve_at(
        control$cens.t.grid, cens.Z[i, , j], time[i], left.limit = TRUE
      )
    }
  }

  obsWeights.event.long <- rep(obsWeights, event.n.time)
  obsWeights.cens.long <- rep(obsWeights, cens.n.time)
  time.event.long <- rep(time, event.n.time)
  time.cens.long <- rep(time, cens.n.time)
  event.event.long <- rep(event, event.n.time)
  event.cens.long <- rep(event, cens.n.time)
  event.t.grid.long <- rep(control$event.t.grid, each = N)
  cens.t.grid.long <- rep(control$cens.t.grid, each = N)

  compute_coef <- function(..., method, event.inclusive = TRUE) {
    .survcomputeCoef(
      ..., method = method, event.inclusive = event.inclusive,
      ipcw.floor = control$ipcw.floor, ipcw.cap = control$ipcw.cap,
      prob.eps = control$logloss.eps,
      optimizer.tol = control$optimizer.tol,
      optimizer.maxit = control$optimizer.maxit
    )
  }

  max_floor_fraction <- 0
  max_cap_fraction <- 0
  min_ipcw_denominator <- Inf
  max_raw_ipcw_weight <- 0
  max_stabilized_ipcw_weight <- 0
  min_ipcw_effective_sample_size <- Inf
  record_diagnostics <- function(weights) {
    diagnostic <- attr(weights, "diagnostics")
    max_floor_fraction <<- max(max_floor_fraction, diagnostic$floor.fraction)
    max_cap_fraction <<- max(max_cap_fraction, diagnostic$cap.fraction)
    weight_diagnostics <- if (!is.null(diagnostic$combined.weights)) {
      diagnostic$combined.weights
    } else {
      diagnostic$weighted.rows
    }
    if (!is.null(weight_diagnostics) && weight_diagnostics$rows > 0) {
      min_ipcw_denominator <<- min(
        min_ipcw_denominator, weight_diagnostics$denominator.minimum
      )
      max_raw_ipcw_weight <<- max(
        max_raw_ipcw_weight, weight_diagnostics$raw.weight.quantiles[["max"]]
      )
      max_stabilized_ipcw_weight <<- max(
        max_stabilized_ipcw_weight,
        weight_diagnostics$stabilized.weight.quantiles[["max"]]
      )
      if (is.finite(weight_diagnostics$effective.sample.size)) {
        min_ipcw_effective_sample_size <<- min(
          min_ipcw_effective_sample_size,
          weight_diagnostics$effective.sample.size
        )
      }
    }
    diagnostic
  }

  initWeightAlg <- get(control$initWeightAlg)
  init.grid <- sort(unique(time))
  obs.cens.vals <- NULL
  obs.event.vals <- NULL
  S.coef <- rep(0, event.k)
  G.coef <- rep(0, cens.k)
  S.grid <- NULL
  G.grid <- NULL
  S.optimizer <- NULL
  G.optimizer <- NULL

  if (control$initWeight == "censoring") {
    initFit <- initWeightAlg(
      time = time, event = 1 - event, X = X, newdata = X,
      new.times = init.grid, obsWeights = obsWeights, id = id
    )
    initial_G <- .row_step_curve_at(
      initFit$pred, init.grid, time, left.limit = TRUE
    )
    obs.cens.vals <- rep(initial_G, event.n.time)
    init_method <- if (metalearner == "logloss") "brier" else metalearner
    initial_coef <- compute_coef(
      time = time.event.long, event = event.event.long,
      t.vals = event.t.grid.long, cens.vals = obs.cens.vals,
      preds = event.Z.long[, !event.errorsInLibrary, drop = FALSE],
      obsWeights = obsWeights.event.long, method = init_method
    )
    S.optimizer <- record_diagnostics(initial_coef)
    S.coef[!event.errorsInLibrary] <- initial_coef
    S.grid <- matrix(drop(event.Z.long %*% S.coef), N, event.n.time)
    obs.event.vals <- rep(
      .row_step_curve_at(S.grid, control$event.t.grid, time), cens.n.time
    )
  } else {
    initFit <- initWeightAlg(
      time = time, event = event, X = X, newdata = X,
      new.times = init.grid, obsWeights = obsWeights, id = id
    )
    initial_S <- .row_step_curve_at(initFit$pred, init.grid, time)
    obs.event.vals <- rep(initial_S, cens.n.time)
  }

  converged <- FALSE
  event.delta <- Inf
  cens.delta <- Inf
  G_t_long <- NULL
  iterations <- 0L

  for (iter in seq_len(control$max.SL.iter)) {
    iterations <- iter
    previous_S.grid <- S.grid
    previous_G.grid <- G.grid

    G.method <- if (metalearner == "logloss") "brier" else metalearner
    censor_coef <- compute_coef(
      time = time.cens.long, event = 1 - event.cens.long,
      t.vals = cens.t.grid.long, cens.vals = obs.event.vals,
      preds = cens.Z.long[, !cens.errorsInLibrary, drop = FALSE],
      obsWeights = obsWeights.cens.long, method = G.method,
      # Censoring arrays represent Pr(C > t); use the matching right-continuous target.
      event.inclusive = TRUE
    )
    G.optimizer <- record_diagnostics(censor_coef)
    G.coef[] <- 0
    G.coef[!cens.errorsInLibrary] <- censor_coef
    G.grid <- matrix(drop(cens.Z.long %*% G.coef), N, cens.n.time)
    obs.cens.vals <- rep(
      .row_step_curve_at(
        G.grid, control$cens.t.grid, time, left.limit = TRUE
      ),
      event.n.time
    )

    if (metalearner == "logloss") {
      G.event.grid <- matrix(NA_real_, N, event.n.time)
      for (i in seq_len(N)) {
        G.event.grid[i, ] <- .step_curve_at(
          control$cens.t.grid, G.grid[i, ], control$event.t.grid
        )
      }
      G_t_long <- as.vector(G.event.grid)
    }

    event_coef <- compute_coef(
      time = time.event.long, event = event.event.long,
      t.vals = event.t.grid.long, cens.vals = obs.cens.vals,
      G_t = if (metalearner == "logloss") G_t_long else NULL,
      preds = event.Z.long[, !event.errorsInLibrary, drop = FALSE],
      obsWeights = obsWeights.event.long, method = metalearner
    )
    S.optimizer <- record_diagnostics(event_coef)
    S.coef[] <- 0
    S.coef[!event.errorsInLibrary] <- event_coef
    S.grid <- matrix(drop(event.Z.long %*% S.coef), N, event.n.time)
    obs.event.vals <- rep(
      .row_step_curve_at(S.grid, control$event.t.grid, time), cens.n.time
    )

    if (!is.null(previous_G.grid) && !is.null(previous_S.grid)) {
      cens.delta <- max(abs(G.grid - previous_G.grid))
      event.delta <- max(abs(S.grid - previous_S.grid))
      converged <- max(cens.delta, event.delta) < control$tol
    }

    if (isTRUE(verbose) || isTRUE(control$traceIter)) {
      message(sprintf(
        "[iter=%d] max grid change: event=%.3e, censoring=%.3e",
        iter, event.delta, cens.delta
      ))
    }
    if (converged) break
  }

  if (!converged) {
    warning("Iterative Super Learner did not converge in ",
            control$max.SL.iter, " iterations.", call. = FALSE)
  } else if (verbose) {
    message("Converged in ", iterations, " iterations.")
  }

  stabilization_fraction <- max(max_floor_fraction, max_cap_fraction)
  if (stabilization_fraction >= control$truncation.warn.fraction &&
      stabilization_fraction > 0) {
    warning(sprintf(
      paste0("IPCW stabilization affected up to %.1f%% of loss rows ",
             "(floor %.1f%%; cap %.1f%%). Consider sensitivity analysis."),
      100 * stabilization_fraction, 100 * max_floor_fraction,
      100 * max_cap_fraction
    ), call. = FALSE)
  }

  calc_risk <- function(preds, time, event, t.grid, cens_T, weights, method,
                        G_t = NULL, event.inclusive = TRUE) {
    indicator <- if (event.inclusive) time <= t.grid else time < t.grid
    stabilized_T <- .stabilize_ipcw(
      cens_T, control$ipcw.floor, control$ipcw.cap
    )
    pseudo_outcome <- 1 - as.numeric(indicator) * event * stabilized_T$weight

    if (method == "brier") {
      return(apply(preds, 2, function(column) {
        mean(weights * (pseudo_outcome - column)^2, na.rm = TRUE)
      }))
    }
    if (method == "entropy") {
      target <- pmax(pmin(pseudo_outcome, 1), 0)
      return(apply(preds, 2, function(column) {
        probability <- pmin(pmax(column, control$logloss.eps),
                            1 - control$logloss.eps)
        mean(-weights * (target * log(probability) +
                           (1 - target) * log(1 - probability)), na.rm = TRUE)
      }))
    }
    if (method == "logloss") {
      if (is.null(G_t)) stop("calc_risk(method='logloss') requires G_t.")
      stabilized_t <- .stabilize_ipcw(
        G_t, control$ipcw.floor, control$ipcw.cap
      )
      fail_weight <- as.numeric(time <= t.grid) * event * stabilized_T$weight
      survive_weight <- as.numeric(time > t.grid) * stabilized_t$weight
      return(apply(preds, 2, function(column) {
        probability <- pmin(pmax(column, control$logloss.eps),
                            1 - control$logloss.eps)
        mean(-weights * (fail_weight * log(1 - probability) +
                           survive_weight * log(probability)), na.rm = TRUE)
      }))
    }
    stop("Unknown method in calc_risk")
  }

  event.cvRisks <- calc_risk(
    event.Z.long, time.event.long, event.event.long, event.t.grid.long,
    obs.cens.vals, obsWeights.event.long, metalearner, G_t_long
  )
  cens.cvRisks <- calc_risk(
    cens.Z.long, time.cens.long, 1 - event.cens.long, cens.t.grid.long,
    obs.event.vals, obsWeights.cens.long, "brier", event.inclusive = TRUE
  )

  list(
    event.coef = S.coef,
    cens.coef = G.coef,
    event.cvRisks = event.cvRisks,
    cens.cvRisks = cens.cvRisks,
    diagnostics = list(
      converged = converged,
      iterations = iterations,
      event.grid.delta = event.delta,
      censoring.grid.delta = cens.delta,
      event.optimizer = S.optimizer,
      censoring.optimizer = G.optimizer,
      ipcw.effective.denominator.threshold = max(
        control$ipcw.floor, 1 / control$ipcw.cap
      ),
      ipcw.minimum.denominator = if (is.finite(min_ipcw_denominator)) {
        min_ipcw_denominator
      } else {
        NA_real_
      },
      ipcw.maximum.raw.weight = max_raw_ipcw_weight,
      ipcw.maximum.stabilized.weight = max_stabilized_ipcw_weight,
      ipcw.minimum.effective.sample.size = if (
        is.finite(min_ipcw_effective_sample_size)
      ) min_ipcw_effective_sample_size else NA_real_,
      ipcw.floor.fraction = max_floor_fraction,
      ipcw.cap.fraction = max_cap_fraction,
      time.weighting = "uniform"
    )
  )
}
