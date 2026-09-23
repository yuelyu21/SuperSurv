# Shared input preparation for optional learner wrappers.
.prepare_wrapper_inputs <- function(time, event, X, newdata, new.times,
                                    obsWeights, learner) {
  if (!is.data.frame(X)) {
    stop("`X` must be a data frame for `", learner, "`.", call. = FALSE)
  }
  if (is.null(newdata)) newdata <- X
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame for `", learner, "`.", call. = FALSE)
  }
  if (!is.numeric(time) || length(time) != nrow(X) ||
      any(!is.finite(time)) || any(time < 0)) {
    stop("`time` must be a finite, non-negative numeric vector with one value per row of `X`.",
         call. = FALSE)
  }
  if (!is.numeric(event) || length(event) != length(time) ||
      any(!event %in% c(0, 1)) || !any(event == 1)) {
    stop("`event` must be a numeric 0/1 vector with at least one observed event.",
         call. = FALSE)
  }
  .validate_data_frame_columns(X, "X")
  .validate_data_frame_columns(newdata, "newdata")
  if (!setequal(names(X), names(newdata))) {
    stop("`newdata` must contain exactly the same feature columns as `X`.",
         call. = FALSE)
  }
  newdata <- newdata[, names(X), drop = FALSE]
  new.times <- .validate_time_grid(new.times, "new.times")

  if (missing(obsWeights) || is.null(obsWeights)) {
    obsWeights <- rep(1, length(time))
  }
  if (!is.numeric(obsWeights) || length(obsWeights) != length(time) ||
      any(!is.finite(obsWeights)) || any(obsWeights < 0) ||
      sum(obsWeights) <= 0) {
    stop("`obsWeights` must be finite, non-negative, and have a positive sum.",
         call. = FALSE)
  }

  list(
    time = as.numeric(time), event = as.numeric(event), X = X,
    newdata = newdata, new.times = new.times,
    obsWeights = as.numeric(obsWeights)
  )
}


.fit_wrapper_matrix_spec <- function(X) {
  factor_levels <- lapply(X[vapply(X, is.factor, logical(1L))], levels)
  matrix <- stats::model.matrix(~ . - 1, data = X)
  list(matrix = matrix, columns = colnames(matrix), factor_levels = factor_levels,
       features = names(X))
}


.predict_wrapper_matrix <- function(newdata, spec, learner) {
  if (!is.null(spec$features)) {
    newdata <- .wrapper_prediction_data(newdata, spec, learner)
  }
  data <- as.data.frame(newdata)
  for (feature in names(spec$factor_levels)) {
    values <- as.character(data[[feature]])
    unexpected <- setdiff(unique(values), spec$factor_levels[[feature]])
    if (length(unexpected)) {
      stop(
        "`newdata` contains unseen level(s) for `", feature, "` in `",
        learner, "`: ", paste(unexpected, collapse = ", "), ".",
        call. = FALSE
      )
    }
    data[[feature]] <- factor(values, levels = spec$factor_levels[[feature]])
  }

  matrix <- stats::model.matrix(~ . - 1, data = data)
  missing_columns <- setdiff(spec$columns, colnames(matrix))
  if (length(missing_columns)) {
    padding <- matrix(0, nrow = nrow(matrix), ncol = length(missing_columns))
    colnames(padding) <- missing_columns
    matrix <- cbind(matrix, padding)
  }
  unexpected_columns <- setdiff(colnames(matrix), spec$columns)
  if (length(unexpected_columns)) {
    stop(
      "`newdata` generated unexpected model-matrix column(s) for `", learner,
      "`: ", paste(unexpected_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  matrix[, spec$columns, drop = FALSE]
}


.finalize_wrapper_survival <- function(prediction, n_observations, times,
                                       learner) {
  prediction <- as.matrix(prediction)
  if (!identical(dim(prediction), c(n_observations, length(times)))) {
    stop(
      "`", learner, "` returned survival predictions with dimensions ",
      paste(dim(prediction), collapse = " x "), "; expected ",
      n_observations, " x ", length(times), ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(prediction))) {
    stop("`", learner, "` returned non-finite survival probabilities.",
         call. = FALSE)
  }
  tolerance <- 1e-8
  if (any(prediction < -tolerance | prediction > 1 + tolerance)) {
    stop("`", learner, "` returned values outside the survival-probability range [0, 1].",
         call. = FALSE)
  }
  prediction <- pmin(pmax(prediction, 0), 1)
  if (ncol(prediction) > 1L) {
    prediction <- t(apply(prediction, 1L, cummin))
  }
  matrix(prediction, nrow = n_observations, ncol = length(times))
}


#' Component-Wise Cox Boosting Learner
#'
#' Fits a weighted component-wise Cox proportional-hazards model using
#' [mboost::glmboost()] and converts its risk score to survival probabilities
#' using SuperSurv's weighted baseline-hazard calibration.
#'
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data frame.
#' @param newdata Covariate data frame used for prediction.
#' @param new.times Times at which survival probabilities are requested.
#' @param obsWeights Optional non-negative observation weights.
#' @param id Currently ignored.
#' @param mstop Number of boosting iterations.
#' @param nu Boosting step size.
#' @param center Whether to center component-wise base learners.
#' @param ties Tied-event approximation for risk-score calibration.
#' @param survival_transform Transformation from calibrated hazard increments
#'   to survival probabilities.
#' @param ... Additional arguments passed to [mboost::glmboost()].
#' @return A list with numeric matrix `pred` and fitted object `fit`.
#' @examples
#' if (requireNamespace("mboost", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:40, ]
#'   X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
#'   fit <- surv.mboost(
#'     dat$duration, dat$event, X, X[1:4, , drop = FALSE],
#'     c(50, 100), mstop = 20
#'   )
#'   dim(fit$pred)
#' }
#' @export
surv.mboost <- function(time, event, X, newdata = NULL, new.times,
                        obsWeights = NULL, id = NULL, mstop = 100L,
                        nu = 0.1, center = FALSE,
                        ties = c("breslow", "efron"),
                        survival_transform = c("exponential", "product_limit"),
                        ...) {
  if (!requireNamespace("mboost", quietly = TRUE)) {
    stop("Learner `surv.mboost` requires the optional package 'mboost'.",
         call. = FALSE)
  }
  input <- .prepare_wrapper_inputs(
    time, event, X, newdata, new.times, obsWeights, "surv.mboost"
  )
  ties <- match.arg(ties)
  survival_transform <- match.arg(survival_transform)
  if (!is.numeric(mstop) || length(mstop) != 1L || !is.finite(mstop) ||
      mstop < 1 || mstop != as.integer(mstop)) {
    stop("`mstop` must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) ||
      nu <= 0 || nu > 1) {
    stop("`nu` must be one number in (0, 1].", call. = FALSE)
  }
  .validate_scalar_logical(center, "center")

  matrix_spec <- .fit_wrapper_matrix_spec(input$X)
  new_matrix <- .predict_wrapper_matrix(input$newdata, matrix_spec, "surv.mboost")
  fit <- mboost::glmboost(
    x = matrix_spec$matrix,
    y = survival::Surv(input$time, input$event),
    family = mboost::CoxPH(),
    weights = input$obsWeights,
    center = center,
    control = mboost::boost_control(mstop = as.integer(mstop), nu = nu),
    ...
  )
  training_score <- as.numeric(stats::predict(
    fit, newdata = matrix_spec$matrix, type = "link"
  ))
  calibration <- .fit_risk_score_calibration(
    input$time, input$event, training_score, input$obsWeights, ties
  )
  new_score <- as.numeric(stats::predict(fit, newdata = new_matrix, type = "link"))
  prediction <- .predict_risk_score_survival(
    calibration, new_score, input$new.times, survival_transform
  )
  prediction <- .finalize_wrapper_survival(
    prediction, nrow(input$newdata), input$new.times, "surv.mboost"
  )

  fit_object <- list(
    object = fit, calibration = calibration, matrix_spec = matrix_spec,
    survival_transform = survival_transform
  )
  class(fit_object) <- "surv.mboost"
  list(pred = prediction, fit = fit_object)
}


#' @noRd
#' @export
predict.surv.mboost <- function(object, newdata, new.times, ...) {
  new.times <- .validate_time_grid(new.times, "new.times")
  matrix <- .predict_wrapper_matrix(newdata, object$matrix_spec, "surv.mboost")
  score <- as.numeric(stats::predict(object$object, newdata = matrix, type = "link"))
  prediction <- .predict_risk_score_survival(
    object$calibration, score, new.times, object$survival_transform
  )
  .finalize_wrapper_survival(prediction, nrow(newdata), new.times, "surv.mboost")
}


.flexsurv_prediction_matrix <- function(object, newdata, new.times,
                                        learner = "surv.flexsurvspline") {
  prediction <- stats::predict(
    object, newdata = newdata, type = "survival", times = new.times,
    conf.int = FALSE
  )
  if (".pred" %in% names(prediction) && is.list(prediction$.pred)) {
    values <- vapply(
      prediction$.pred,
      function(value) as.numeric(value$.pred_survival),
      numeric(length(new.times))
    )
    prediction <- t(matrix(values, nrow = length(new.times), ncol = nrow(newdata)))
  } else if (".pred_survival" %in% names(prediction)) {
    prediction <- matrix(
      prediction$.pred_survival,
      nrow = nrow(newdata), ncol = length(new.times), byrow = TRUE
    )
  } else {
    stop("`", learner, "` could not extract survival probabilities from 'flexsurv'.",
         call. = FALSE)
  }
  .finalize_wrapper_survival(
    prediction, nrow(newdata), new.times, learner
  )
}


#' Flexible Parametric Spline Survival Learner
#'
#' Fits a weighted Royston-Parmar flexible parametric survival model using
#' [flexsurv::flexsurvspline()].
#'
#' @inheritParams surv.mboost
#' @param k Number of internal spline knots.
#' @param scale Scale on which the flexible model is defined.
#' @param ... Additional arguments passed to [flexsurv::flexsurvspline()].
#' @return A list with numeric matrix `pred` and fitted object `fit`.
#' @examples
#' if (requireNamespace("flexsurv", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:40, ]
#'   X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
#'   fit <- surv.flexsurvspline(
#'     dat$duration, dat$event, X, X[1:4, , drop = FALSE],
#'     c(50, 100), k = 1
#'   )
#'   dim(fit$pred)
#' }
#' @export
surv.flexsurvspline <- function(time, event, X, newdata = NULL, new.times,
                                obsWeights = NULL, id = NULL, k = 1L,
                                scale = "hazard", ...) {
  if (!requireNamespace("flexsurv", quietly = TRUE)) {
    stop("Learner `surv.flexsurvspline` requires the optional package 'flexsurv'.",
         call. = FALSE)
  }
  input <- .prepare_wrapper_inputs(
    time, event, X, newdata, new.times, obsWeights, "surv.flexsurvspline"
  )
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 0 || k != as.integer(k)) {
    stop("`k` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.character(scale) || length(scale) != 1L || is.na(scale)) {
    stop("`scale` must be one character string.", call. = FALSE)
  }

  training_data <- data.frame(
    .supersurv_time = input$time,
    .supersurv_event = input$event,
    input$X,
    check.names = FALSE
  )
  formula <- stats::as.formula(
    "survival::Surv(.supersurv_time, .supersurv_event) ~ ."
  )
  fit <- flexsurv::flexsurvspline(
    formula = formula, data = training_data, weights = input$obsWeights,
    k = as.integer(k), scale = scale, ...
  )
  prediction <- .flexsurv_prediction_matrix(
    fit, input$newdata, input$new.times
  )

  fit_object <- list(object = fit, features = names(input$X))
  class(fit_object) <- "surv.flexsurvspline"
  list(pred = prediction, fit = fit_object)
}


#' @noRd
#' @export
predict.surv.flexsurvspline <- function(object, newdata, new.times, ...) {
  if (!is.data.frame(newdata) || !setequal(names(newdata), object$features)) {
    stop("`newdata` must contain exactly the features used by `surv.flexsurvspline`.",
         call. = FALSE)
  }
  newdata <- newdata[, object$features, drop = FALSE]
  new.times <- .validate_time_grid(new.times, "new.times")
  .flexsurv_prediction_matrix(object$object, newdata, new.times)
}


#' Generalized Random Survival Forest Learner
#'
#' Fits an honest generalized random survival forest using
#' [grf::survival_forest()] and returns conditional survival curves.
#'
#' @inheritParams surv.mboost
#' @param num.trees Number of trees.
#' @param mtry Number of candidate variables considered at each split.
#' @param min.node.size Minimum terminal-node size.
#' @param honesty Whether to use honest sample splitting.
#' @param prediction.type Either `"Kaplan-Meier"` or `"Nelson-Aalen"`.
#' @param seed Integer random seed passed to `grf`.
#' @param ... Additional arguments passed to [grf::survival_forest()].
#' @return A list with numeric matrix `pred` and fitted object `fit`.
#' @examples
#' if (requireNamespace("grf", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:40, ]
#'   X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
#'   fit <- surv.grf(
#'     dat$duration, dat$event, X, X[1:4, , drop = FALSE],
#'     c(50, 100), num.trees = 50, seed = 1
#'   )
#'   dim(fit$pred)
#' }
#' @export
surv.grf <- function(time, event, X, newdata = NULL, new.times,
                     obsWeights = NULL, id = NULL, num.trees = 1000L,
                     mtry = NULL, min.node.size = 15L, honesty = TRUE,
                     prediction.type = c("Kaplan-Meier", "Nelson-Aalen"),
                     seed = 1L, ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Learner `surv.grf` requires the optional package 'grf'.",
         call. = FALSE)
  }
  input <- .prepare_wrapper_inputs(
    time, event, X, newdata, new.times, obsWeights, "surv.grf"
  )
  prediction.type <- match.arg(prediction.type)
  .validate_scalar_logical(honesty, "honesty")
  for (argument in c("num.trees", "min.node.size", "seed")) {
    value <- get(argument)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      stop("`", argument, "` must be one positive integer.", call. = FALSE)
    }
  }

  matrix_spec <- .fit_wrapper_matrix_spec(input$X)
  new_matrix <- .predict_wrapper_matrix(input$newdata, matrix_spec, "surv.grf")
  arguments <- list(
    X = matrix_spec$matrix, Y = input$time, D = input$event,
    num.trees = as.integer(num.trees), sample.weights = input$obsWeights,
    min.node.size = as.integer(min.node.size), honesty = honesty,
    prediction.type = prediction.type, seed = as.integer(seed)
  )
  if (!is.null(mtry)) arguments$mtry <- mtry
  arguments <- c(arguments, list(...))
  fit <- do.call(grf::survival_forest, arguments)
  prediction <- stats::predict(
    fit, newdata = new_matrix, failure.times = input$new.times,
    prediction.times = "curve", prediction.type = prediction.type
  )$predictions
  prediction <- .finalize_wrapper_survival(
    prediction, nrow(input$newdata), input$new.times, "surv.grf"
  )

  fit_object <- list(
    object = fit, matrix_spec = matrix_spec, prediction.type = prediction.type
  )
  class(fit_object) <- "surv.grf"
  list(pred = prediction, fit = fit_object)
}


#' @noRd
#' @export
predict.surv.grf <- function(object, newdata, new.times, ...) {
  new.times <- .validate_time_grid(new.times, "new.times")
  matrix <- .predict_wrapper_matrix(newdata, object$matrix_spec, "surv.grf")
  prediction <- stats::predict(
    object$object, newdata = matrix, failure.times = new.times,
    prediction.times = "curve", prediction.type = object$prediction.type
  )$predictions
  .finalize_wrapper_survival(prediction, nrow(newdata), new.times, "surv.grf")
}


.require_pycox_backend <- function(learner) {
  if (!requireNamespace("survivalmodels", quietly = TRUE)) {
    stop("Learner `", learner,
         "` requires the optional R package 'survivalmodels'.", call. = FALSE)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Learner `", learner,
         "` requires the optional R package 'reticulate'.", call. = FALSE)
  }
  modules <- c("torch", "torchtuples", "pycox")
  available <- vapply(modules, reticulate::py_module_available, logical(1L))
  if (!all(available)) {
    stop(
      "Learner `", learner, "` could not import the required Python module(s) ",
      paste(sprintf("'%s'", modules[!available]), collapse = ", "),
      " from the Python environment used by 'reticulate'. Install or repair ",
      "the module(s) before fitting; SuperSurv does not install Python software automatically.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


.validate_deep_wrapper_arguments <- function(input, learner, num_nodes,
                                             activation, batch_norm, dropout,
                                             epochs, batch_size, verbose,
                                             seed) {
  if (length(unique(input$obsWeights)) > 1L) {
    stop(
      "`", learner, "` does not support nonuniform `obsWeights` because its ",
      "'survivalmodels'/'pycox' backend does not expose case weights. Use ",
      "uniform weights or choose a learner with native weight support.",
      call. = FALSE
    )
  }
  if (!is.numeric(num_nodes) || !length(num_nodes) ||
      any(!is.finite(num_nodes)) || any(num_nodes < 1) ||
      any(num_nodes != as.integer(num_nodes))) {
    stop("`num_nodes` must contain positive integers.", call. = FALSE)
  }
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !nzchar(activation)) {
    stop("`activation` must be one non-empty character string.", call. = FALSE)
  }
  .validate_scalar_logical(batch_norm, "batch_norm")
  .validate_scalar_logical(verbose, "verbose")
  if (!is.null(dropout) &&
      (!is.numeric(dropout) || length(dropout) != 1L ||
       !is.finite(dropout) || dropout < 0 || dropout >= 1)) {
    stop("`dropout` must be NULL or one number in [0, 1).", call. = FALSE)
  }
  for (argument in c("epochs", "batch_size", "seed")) {
    value <- get(argument)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      stop("`", argument, "` must be one positive integer.", call. = FALSE)
    }
  }
  invisible(TRUE)
}


.set_pycox_seed <- function(seed) {
  code <- paste0(
    "import random, numpy as np, torch\n",
    "random.seed(", as.integer(seed), ")\n",
    "np.random.seed(", as.integer(seed), ")\n",
    "torch.manual_seed(", as.integer(seed), ")\n",
    "if torch.cuda.is_available(): torch.cuda.manual_seed_all(",
    as.integer(seed), ")"
  )
  reticulate::py_run_string(code)
  invisible(NULL)
}


.deep_survival_matrix <- function(object, newdata, new.times, learner,
                                  batch_size, ...) {
  if (reticulate::py_is_null_xptr(object$model)) {
    stop("`", learner, "` no longer has a live Python fitted model. ",
         "Refit in the active R/Python session; plain saveRDS()/readRDS() ",
         "does not preserve Python model objects.", call. = FALSE)
  }
  native <- stats::predict(
    object, newdata = as.data.frame(newdata), type = "survival",
    batch_size = as.integer(batch_size), ...
  )
  native <- as.matrix(native)
  native_times <- suppressWarnings(as.numeric(colnames(native)))
  if (ncol(native) < 1L || length(native_times) != ncol(native) ||
      any(!is.finite(native_times))) {
    stop("`", learner,
         "` could not recover the native prediction-time grid from 'survivalmodels'.",
         call. = FALSE)
  }
  ordering <- order(native_times)
  native_times <- native_times[ordering]
  native <- native[, ordering, drop = FALSE]
  keep <- !duplicated(native_times, fromLast = TRUE)
  native_times <- native_times[keep]
  native <- native[, keep, drop = FALSE]

  indices <- findInterval(new.times, native_times)
  prediction <- matrix(1, nrow = nrow(native), ncol = length(new.times))
  after_first <- indices > 0L
  if (any(after_first)) {
    prediction[, after_first] <- native[, indices[after_first], drop = FALSE]
  }
  .finalize_wrapper_survival(prediction, nrow(newdata), new.times, learner)
}


#' Experimental DeepSurv Learner
#'
#' Fits a Cox partial-likelihood neural network through
#' [survivalmodels::deepsurv()]. This optional adapter requires a Python
#' environment containing `torch`, `torchtuples`, and `pycox` and currently
#' supports only uniform observation weights.
#'
#' @inheritParams surv.mboost
#' @param num_nodes Positive integers giving hidden-layer sizes.
#' @param activation Neural-network activation name.
#' @param batch_norm Whether to use batch normalization.
#' @param dropout Optional dropout probability in `[0, 1)`.
#' @param epochs Number of training epochs.
#' @param batch_size Training and prediction batch size.
#' @param device Optional device passed to `survivalmodels`.
#' @param verbose Whether the Python backend should print training progress.
#' @param seed Positive integer used for Python, NumPy, and torch random-number
#'   generators.
#' @param ... Additional arguments passed to [survivalmodels::deepsurv()].
#' @details Fitted Python objects can be reused in the active R/Python session.
#'   Plain `saveRDS()` is not a portable persistence format for these objects.
#'   Native survival curves are mapped to requested times as right-continuous
#'   steps; see [surv.coxtime()] for the boundary convention.
#' @return A list with numeric matrix `pred` and fitted object `fit`.
#' @examples
#' if (interactive() && requireNamespace("survivalmodels", quietly = TRUE) &&
#'     requireNamespace("reticulate", quietly = TRUE) &&
#'     reticulate::py_module_available("pycox")) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:60, ]
#'   X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
#'   fit <- surv.deepsurv(
#'     dat$duration, dat$event, X, X[1:4, , drop = FALSE],
#'     c(50, 100), epochs = 2, seed = 1
#'   )
#'   dim(fit$pred)
#' }
#' @export
surv.deepsurv <- function(time, event, X, newdata = NULL, new.times,
                          obsWeights = NULL, id = NULL,
                          num_nodes = c(32L, 32L), activation = "relu",
                          batch_norm = TRUE, dropout = NULL,
                          epochs = 100L, batch_size = 128L,
                          device = NULL, verbose = FALSE, seed = 1L, ...) {
  .require_pycox_backend("surv.deepsurv")
  input <- .prepare_wrapper_inputs(
    time, event, X, newdata, new.times, obsWeights, "surv.deepsurv"
  )
  .validate_deep_wrapper_arguments(
    input, "surv.deepsurv", num_nodes, activation, batch_norm, dropout,
    epochs, batch_size, verbose, seed
  )
  matrix_spec <- .fit_wrapper_matrix_spec(input$X)
  new_matrix <- .predict_wrapper_matrix(input$newdata, matrix_spec, "surv.deepsurv")
  .set_pycox_seed(seed)
  fit <- survivalmodels::deepsurv(
    x = as.data.frame(matrix_spec$matrix),
    y = survival::Surv(input$time, input$event),
    num_nodes = as.integer(num_nodes), activation = activation,
    batch_norm = batch_norm, dropout = dropout, device = device,
    epochs = as.integer(epochs), batch_size = as.integer(batch_size),
    verbose = verbose, ...
  )
  prediction <- .deep_survival_matrix(
    fit, new_matrix, input$new.times, "surv.deepsurv", batch_size
  )

  fit_object <- list(
    object = fit, matrix_spec = matrix_spec, batch_size = as.integer(batch_size)
  )
  class(fit_object) <- "surv.deepsurv"
  list(pred = prediction, fit = fit_object)
}


#' @noRd
#' @export
predict.surv.deepsurv <- function(object, newdata, new.times,
                                  batch_size = object$batch_size, ...) {
  .require_pycox_backend("surv.deepsurv")
  new.times <- .validate_time_grid(new.times, "new.times")
  matrix <- .predict_wrapper_matrix(newdata, object$matrix_spec, "surv.deepsurv")
  .deep_survival_matrix(
    object$object, matrix, new.times, "surv.deepsurv", batch_size, ...
  )
}


#' Experimental DeepHit Learner
#'
#' Fits a single-event discrete-time DeepHit neural network through
#' [survivalmodels::deephit()]. This optional adapter requires a Python
#' environment containing `torch`, `torchtuples`, and `pycox` and currently
#' supports only uniform observation weights.
#'
#' @inheritParams surv.deepsurv
#' @param cuts Number of discrete time intervals, or a vector of cut points
#'   supplied through `cutpoints`.
#' @param cutpoints Optional numeric vector of cut points.
#' @param scheme Cut-point scheme used when `cutpoints` is `NULL`.
#' @param mod_alpha Weight placed on the likelihood component of the DeepHit
#'   objective.
#' @param sigma Ranking-loss bandwidth.
#' @param ... Additional arguments passed to [survivalmodels::deephit()].
#' @details Fitted Python objects can be reused in the active R/Python session.
#'   Plain `saveRDS()` is not a portable persistence format for these objects.
#'   Native survival curves are mapped to requested times as right-continuous
#'   steps; see [surv.coxtime()] for the boundary convention.
#' @return A list with numeric matrix `pred` and fitted object `fit`.
#' @examples
#' if (interactive() && requireNamespace("survivalmodels", quietly = TRUE) &&
#'     requireNamespace("reticulate", quietly = TRUE) &&
#'     reticulate::py_module_available("pycox")) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:60, ]
#'   X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
#'   fit <- surv.deephit(
#'     dat$duration, dat$event, X, X[1:4, , drop = FALSE],
#'     c(50, 100), cuts = 10, epochs = 2, seed = 1
#'   )
#'   dim(fit$pred)
#' }
#' @export
surv.deephit <- function(time, event, X, newdata = NULL, new.times,
                         obsWeights = NULL, id = NULL,
                         num_nodes = c(32L, 32L), activation = "relu",
                         batch_norm = TRUE, dropout = NULL,
                         epochs = 100L, batch_size = 128L,
                         device = NULL, verbose = FALSE, seed = 1L,
                         cuts = 20L, cutpoints = NULL,
                         scheme = c("equidistant", "quantiles"),
                         mod_alpha = 0.2, sigma = 0.1, ...) {
  .require_pycox_backend("surv.deephit")
  input <- .prepare_wrapper_inputs(
    time, event, X, newdata, new.times, obsWeights, "surv.deephit"
  )
  .validate_deep_wrapper_arguments(
    input, "surv.deephit", num_nodes, activation, batch_norm, dropout,
    epochs, batch_size, verbose, seed
  )
  scheme <- match.arg(scheme)
  if (is.null(cutpoints)) {
    if (!is.numeric(cuts) || length(cuts) != 1L || !is.finite(cuts) ||
        cuts < 2 || cuts != as.integer(cuts)) {
      stop("`cuts` must be one integer of at least 2.", call. = FALSE)
    }
  } else if (!is.numeric(cutpoints) || length(cutpoints) < 2L ||
             any(!is.finite(cutpoints)) || any(diff(cutpoints) <= 0)) {
    stop("`cutpoints` must be a strictly increasing finite numeric vector.",
         call. = FALSE)
  }
  if (!is.numeric(mod_alpha) || length(mod_alpha) != 1L ||
      !is.finite(mod_alpha) || mod_alpha < 0 || mod_alpha > 1) {
    stop("`mod_alpha` must be one number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L ||
      !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be one positive number.", call. = FALSE)
  }

  matrix_spec <- .fit_wrapper_matrix_spec(input$X)
  new_matrix <- .predict_wrapper_matrix(input$newdata, matrix_spec, "surv.deephit")
  .set_pycox_seed(seed)
  fit <- survivalmodels::deephit(
    x = as.data.frame(matrix_spec$matrix),
    y = survival::Surv(input$time, input$event),
    num_nodes = as.integer(num_nodes), activation = activation,
    batch_norm = batch_norm, dropout = dropout, device = device,
    epochs = as.integer(epochs), batch_size = as.integer(batch_size),
    verbose = verbose, cuts = as.integer(cuts), cutpoints = cutpoints,
    scheme = scheme, mod_alpha = mod_alpha, sigma = sigma, ...
  )
  prediction <- .deep_survival_matrix(
    fit, new_matrix, input$new.times, "surv.deephit", batch_size
  )

  fit_object <- list(
    object = fit, matrix_spec = matrix_spec, batch_size = as.integer(batch_size)
  )
  class(fit_object) <- "surv.deephit"
  list(pred = prediction, fit = fit_object)
}


#' @noRd
#' @export
predict.surv.deephit <- function(object, newdata, new.times,
                                 batch_size = object$batch_size, ...) {
  .require_pycox_backend("surv.deephit")
  new.times <- .validate_time_grid(new.times, "new.times")
  matrix <- .predict_wrapper_matrix(newdata, object$matrix_spec, "surv.deephit")
  .deep_survival_matrix(
    object$object, matrix, new.times, "surv.deephit", batch_size, ...
  )
}


#' Experimental Cox-Time Neural Survival Learner
#'
#' Fits a neural Cox-Time model allowing time-dependent covariate effects through
#' [survivalmodels::coxtime()]. This optional adapter requires Python `torch`,
#' `torchtuples`, and `pycox` and supports only uniform observation weights.
#'
#' @inheritParams surv.deepsurv
#' @param standardize_time Whether to standardize outcome times using the
#'   training-only Cox-Time label transformation. Predictions are returned on
#'   the original time scale.
#' @param ... Additional named fitting arguments passed to
#'   [survivalmodels::coxtime()], such as `frac` for an internal training-fold
#'   validation split or `early_stopping`. Outcomes and features are managed
#'   by the adapter.
#' @details Native survival curves are evaluated as right-continuous steps on
#'   the backend's time grid, using survival one before the first grid point
#'   and the final available value beyond the last point. This is not a claim
#'   of reliable extrapolation beyond observed follow-up. Python-backed fitted
#'   objects are intended for reuse within the active R/Python session; plain
#'   `saveRDS()` is not a portable persistence format for Python objects.
#' @return A list with numeric survival matrix `pred` and fitted object `fit`.
#' @examples
#' if (interactive() && requireNamespace("survivalmodels", quietly = TRUE) &&
#'     requireNamespace("reticulate", quietly = TRUE) &&
#'     reticulate::py_module_available("pycox")) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:60, ]
#'   X <- dat[, "x1", drop = FALSE]
#'   fit <- surv.coxtime(dat$duration, dat$event, X,
#'                       X[1:3, , drop = FALSE], c(50, 100),
#'                       epochs = 2, batch_norm = FALSE, device = "cpu")
#'   dim(fit$pred)
#' }
#' @export
surv.coxtime <- function(time, event, X, newdata = NULL, new.times,
                         obsWeights = NULL, id = NULL,
                         num_nodes = c(32L, 32L), activation = "relu",
                         batch_norm = TRUE, dropout = NULL,
                         epochs = 100L, batch_size = 128L,
                         device = NULL, verbose = FALSE, seed = 1L,
                         standardize_time = TRUE, ...) {
  learner <- "surv.coxtime"
  input <- .prepare_wrapper_inputs(time, event, X, newdata, new.times,
                                   obsWeights, learner)
  .validate_deep_wrapper_arguments(input, learner, num_nodes, activation,
                                    batch_norm, dropout, epochs, batch_size,
                                    verbose, seed)
  .validate_scalar_logical(standardize_time, "standardize_time")
  dots <- list(...)
  .check_wrapper_dots(dots, c("formula", "data", "x", "y", "reverse",
                             "time_variable", "status_variable"), learner)
  .require_pycox_backend(learner)
  matrix_spec <- .fit_wrapper_matrix_spec(input$X)
  new_matrix <- .predict_wrapper_matrix(input$newdata, matrix_spec, learner)
  .set_pycox_seed(seed)
  fit <- do.call(survivalmodels::coxtime, c(list(
    x = as.data.frame(matrix_spec$matrix),
    y = survival::Surv(input$time, input$event),
    num_nodes = as.integer(num_nodes), activation = activation,
    batch_norm = batch_norm, dropout = dropout, device = device,
    epochs = as.integer(epochs), batch_size = as.integer(batch_size),
    verbose = verbose, standardize_time = standardize_time
  ), dots))
  prediction <- .deep_survival_matrix(fit, new_matrix, input$new.times,
                                      learner, batch_size)
  fit_object <- list(object = fit, matrix_spec = matrix_spec,
                     batch_size = as.integer(batch_size))
  class(fit_object) <- learner
  list(pred = prediction, fit = fit_object)
}

#' @noRd
#' @export
predict.surv.coxtime <- function(object, newdata, new.times,
                                batch_size = object$batch_size, ...) {
  .require_pycox_backend("surv.coxtime")
  new.times <- .validate_time_grid(new.times, "new.times")
  matrix <- .predict_wrapper_matrix(newdata, object$matrix_spec, "surv.coxtime")
  .deep_survival_matrix(object$object, matrix, new.times, "surv.coxtime",
                        batch_size, ...)
}
