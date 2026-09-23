#' Create a Tuning Grid of Survival Learners
#'
#' Dynamically generates custom wrapper functions for a specified base learner
#' across a grid of hyperparameters.
#'
#' @param base_learner Character string of the base learner function name (e.g., "surv.gbm").
#' @param grid_params List of numeric/character vectors containing hyperparameter values.
#' @return A character vector of class \code{"SuperSurv_grid"} containing the
#'   newly generated function names.
#' @details Generated functions are assigned to the calling environment. The
#'   base learner is resolved there (with package functions as a fallback) when
#'   the grid is created and retained by each generated function. Local base
#'   learners and grids can therefore be used inside a function that calls
#'   \code{SuperSurv()}, without modifying the global environment.
#' @export
#' @keywords internal
create_grid <- function(base_learner, grid_params) {
  if (!is.character(base_learner) || length(base_learner) != 1L ||
      is.na(base_learner) || !nzchar(base_learner)) {
    stop("`base_learner` must be one non-empty function name.", call. = FALSE)
  }
  base_fun <- .find_library_function(base_learner, parent.frame())
  if (is.null(base_fun)) {
    stop("`base_learner` contains an unknown function name: ", base_learner,
         ". Define the function before calling `create_grid()`.", call. = FALSE)
  }
  param_grid <- expand.grid(grid_params, stringsAsFactors = FALSE)

  generated_learners <- character(nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {

    specific_params <- as.list(param_grid[i, , drop = FALSE])

    # Create descriptive learner name
    param_string <- paste(
      paste0(names(specific_params), specific_params),
      collapse = "_"
    )

    learner_name <- paste0(base_learner, "_", param_string)

    generated_learners[i] <- learner_name

    new_fun <- local({

      .base_fun <- base_fun
      .specific_params <- specific_params

      function(time, event, X, newdata, new.times, obsWeights = NULL, id = NULL, ...) {

        args_to_pass <- list(
          time = time,
          event = event,
          X = X,
          newdata = newdata,
          new.times = new.times,
          obsWeights = obsWeights,
          id = id
        )

        call_args <- c(args_to_pass, .specific_params, list(...))

        do.call(.base_fun, call_args)
      }

    })

    assign(learner_name, new_fun, envir = parent.frame())
  }

  class(generated_learners) <- c("SuperSurv_grid", class(generated_learners))
  return(generated_learners)
}






#' Fit a baseline cumulative-hazard calibration for a fixed risk score
#' @noRd
.fit_risk_score_calibration <- function(time, event, risk_score,
                                        obsWeights = NULL,
                                        ties = c("breslow", "efron")) {
  ties <- match.arg(ties)
  if (is.null(obsWeights)) obsWeights <- rep(1, length(time))
  if (!length(time) || length(event) != length(time) ||
      length(risk_score) != length(time) || length(obsWeights) != length(time)) {
    stop("Calibration inputs must have the same positive length.", call. = FALSE)
  }
  if (any(!is.finite(time)) || any(!is.finite(risk_score)) ||
      any(!is.finite(obsWeights)) || any(obsWeights < 0)) {
    stop("Calibration times, risk scores, and weights must be finite; weights must be non-negative.",
         call. = FALSE)
  }
  if (any(!event %in% c(0, 1))) {
    stop("Calibration event indicators must be zero or one.", call. = FALSE)
  }
  if (!any(event == 1)) {
    return(list(
      event.times = numeric(), hazard.increment = numeric(),
      cumulative.hazard = numeric(), center = stats::weighted.mean(
        risk_score, obsWeights
      ), ties = ties
    ))
  }

  center <- stats::weighted.mean(risk_score, obsWeights)
  centered_score <- risk_score - center
  event_times <- sort(unique(time[event == 1]))

  if (ties == "breslow") {
    hazard_increment <- vapply(event_times, function(event_time) {
      weighted_events <- sum(obsWeights[time == event_time & event == 1])
      risk_set <- time >= event_time
      denominator <- sum(obsWeights[risk_set] * exp(centered_score[risk_set]))
      weighted_events / denominator
    }, numeric(1L))
    cumulative_hazard <- cumsum(hazard_increment)
  } else {
    calibration_fit <- survival::coxph(
      survival::Surv(time, event) ~ offset(centered_score),
      weights = obsWeights,
      ties = "efron",
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
    baseline <- survival::basehaz(calibration_fit, centered = FALSE)
    cumulative_hazard <- vapply(event_times, function(event_time) {
      values <- baseline$hazard[baseline$time <= event_time]
      if (length(values)) utils::tail(values, 1L) else 0
    }, numeric(1L))
    hazard_increment <- diff(c(0, cumulative_hazard))
  }

  list(
    event.times = event_times,
    hazard.increment = pmax(hazard_increment, 0),
    cumulative.hazard = cummax(pmax(cumulative_hazard, 0)),
    center = center,
    ties = ties
  )
}


#' Convert calibrated proportional-hazards risk scores to survival curves
#' @noRd
.predict_risk_score_survival <- function(calibration, risk_score, new.times,
                                         survival_transform = c(
                                           "exponential", "product_limit"
                                         )) {
  survival_transform <- match.arg(survival_transform)
  centered_score <- pmin(risk_score - calibration$center, 700)
  relative_risk <- exp(centered_score)

  if (!length(calibration$event.times)) {
    return(matrix(1, nrow = length(risk_score), ncol = length(new.times)))
  }

  if (survival_transform == "exponential") {
    cumulative_hazard <- vapply(new.times, function(eval_time) {
      index <- which(calibration$event.times <= eval_time)
      if (length(index)) calibration$cumulative.hazard[max(index)] else 0
    }, numeric(1L))
    prediction <- outer(
      relative_risk, cumulative_hazard,
      function(risk, hazard) exp(-risk * hazard)
    )
  } else {
    prediction <- vapply(new.times, function(eval_time) {
      index <- which(calibration$event.times <= eval_time)
      if (!length(index)) return(rep(1, length(relative_risk)))
      vapply(relative_risk, function(risk) {
        factors <- pmax(1 - risk * calibration$hazard.increment[index], 0)
        prod(factors)
      }, numeric(1L))
    }, numeric(length(relative_risk)))
    prediction <- matrix(
      prediction, nrow = length(relative_risk), ncol = length(new.times)
    )
  }

  prediction <- pmin(pmax(prediction, 0), 1)
  if (ncol(prediction) > 1L) prediction <- t(apply(prediction, 1, cummin))
  prediction
}


#' Calculate the Breslow estimator of the cumulative baseline hazard
#' @noRd
safe_breslow_step <- function(time, event, risk_score, new.times,
                              obsWeights = NULL) {
  calibration <- .fit_risk_score_calibration(
    time, event, risk_score, obsWeights, ties = "breslow"
  )
  if (!length(calibration$event.times)) return(rep(0, length(new.times)))
  vapply(new.times, function(eval_time) {
    index <- which(calibration$event.times <= eval_time)
    if (length(index)) calibration$cumulative.hazard[max(index)] else 0
  }, numeric(1L))
}







#' Universal Risk Wrapper for SHAP
#'
#' Ensures "High Value = High Risk" orientation for ALL base learners and the
#' SuperSurv ensemble, allowing unified SHAP value calculations.
#'
#' @param object A fitted model object (e.g., from a single wrapper or a SuperSurv ensemble).
#' @param newdata A data.frame of new covariates to predict on.
#'
#' @return A numeric vector of risk scores of the same length as the number
#'   of rows in \code{newdata}, where higher values consistently indicate a higher
#'   risk of the event.
#' @keywords internal
#' @noRd
get_risk_universal <- function(object, newdata) {

  # --- 1. The Cox Family (High = Bad) ---
  if (inherits(object, "surv.coxph") || inherits(object, "coxph")) {
    # Extract native model to avoid predict.surv.coxph collision
    mod <- if(inherits(object, "surv.coxph")) object$object else object
    return(predict(mod, newdata = newdata, type = "lp"))
  }
  # --- REPLACE THIS BLOCK IN get_risk_universal ---
  else if (inherits(object, "surv.glmnet")) {
    mod <- if(inherits(object, "surv.glmnet")) object$object else object

    # glmnet REQUIRES a matrix, so we convert the SHAP background data
    newdata_df <- as.data.frame(newdata)
    newdata_mat <- stats::model.matrix(~ . - 1, data = newdata_df)

    # Predict the Linear Predictor (Risk)
    return(as.vector(predict(mod, newx = newdata_mat, s = "lambda.min", type = "link")))
  }
  else if (inherits(object, "surv.coxboost")) {
    mod <- object$object
    return(as.vector(predict(mod, newdata = as.matrix(newdata), type = "lp")))
  }

  # --- 2. The GAM ---
  else if (inherits(object, "surv.gam") || inherits(object, "gam")) {
    mod <- if(inherits(object, "surv.gam")) object$object else object
    return(as.vector(predict(mod, newdata = newdata, type = "link")))
  }

  # --- 3. The Machine Family (High = Bad) ---
  else if (inherits(object, "surv.gbm") || inherits(object, "gbm")) {
    mod <- if(inherits(object, "surv.gbm")) object$object else object
    best_iter <- if(inherits(object, "surv.gbm")) object$best.iter else mod$n.trees
    return(predict(mod, newdata = newdata, n.trees = best_iter, type = "link"))
  }
  else if (inherits(object, "surv.xgboost") || inherits(object, "xgb.Booster")) {
    mod <- if(inherits(object, "surv.xgboost")) object$object else object
    if(is.data.frame(newdata)) newdata <- as.matrix(newdata)
    return(predict(mod, newdata = newdata, outputmargin = TRUE))
  }
  else if (inherits(object, "surv.svm")) {
    mod <- object$object
    return(predict(mod, newdata = newdata)$predicted)
  }

  # --- 4. The Tree Family (High = Bad via Conversion) ---
  else if (inherits(object, "surv.rfsrc") || inherits(object, "rfsrc")) {
    mod <- if(inherits(object, "surv.rfsrc")) object$object else object
    return(predict(mod, newdata = newdata)$predicted)
  }
  else if (inherits(object, "surv.ranger") || inherits(object, "ranger")) {
    mod <- if(inherits(object, "surv.ranger")) object$object else object
    chf <- predict(mod, data = newdata)$chf
    mid_idx <- floor(ncol(chf)/2)
    return(chf[, mid_idx])
  }
  else if (inherits(object, "surv.bart")) {
    mod <- object$object
    surv_mat <- mod$surv.test.mean
    return(1 - rowMeans(surv_mat))
  }
  else if (inherits(object, "surv.aorsf") || inherits(object, "aorsf")) {
    mod <- if(inherits(object, "surv.aorsf")) object$object else object

    # 1. Get the Median Time Horizon (from training data)
    t_median <- median(mod$data[[1]], na.rm = TRUE)

    # 2. Predict Risk Directly (Fastest & Safest)
    risk <- predict(
      mod,
      new_data = newdata,
      pred_horizon = t_median,
      pred_type = "risk"
    )

    return(as.numeric(risk))
  }
  # Modify this line in get_risk_universal:
  else if (inherits(object, "surv.glmnet") || inherits(object, "surv.ridge")) {

    mod <- if(inherits(object, "surv.glmnet") || inherits(object, "surv.ridge")) object$object else object

    newdata_df <- as.data.frame(newdata)
    newdata_mat <- stats::model.matrix(~ . - 1, data = newdata_df)

    return(as.vector(predict(mod, newx = newdata_mat, s = "lambda.min", type = "link")))
  }

  # --- 5. The Parametric Family (High = GOOD -> FLIP SIGN!) ---
  else if (inherits(object, "surv.parametric") ||
           inherits(object, "survreg") ||
           inherits(object, "surv.weibull") ||
           inherits(object, "surv.exponential") ||
           inherits(object, "surv.lognormal") ||
           inherits(object, "surv.loglogistic")) {

    mod <- object$object
    # We use the linear predictor as the risk score for SHAP
    return(as.vector(predict(mod, newdata = as.data.frame(newdata), type = "linear")))
  }
  else if (inherits(object, "surv.rpart")) {
    mod <- object$object

    # rpart for survival returns the expected event rate.
    # We take the log (with a safety clamp) to get the linear predictor (risk score) for SHAP.
    rate_new <- predict(mod, newdata = as.data.frame(newdata))
    return(as.vector(log(pmax(rate_new, 1e-10))))
  }

  else {
    stop(paste("Unknown class in SHAP wrapper:", class(object)[1]))
  }
}
