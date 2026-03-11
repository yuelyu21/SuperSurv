#' Create a Tuning Grid of Survival Learners
#'
#' Dynamically generates custom wrapper functions for a specified base learner
#' across a grid of hyperparameters.
#'
#' @param base_learner Character string of the base learner function name (e.g., "surv.gbm").
#' @param grid_params List of numeric/character vectors containing hyperparameter values.
#' @return A character vector of the newly generated function names.
#' @export
#' @keywords internal
create_grid <- function(base_learner, grid_params) {

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

      .base_learner <- base_learner
      .specific_params <- specific_params

      function(time, event, X, newdata, new.times, obsWeights = NULL, id = NULL, ...) {

        base_fun <- get(.base_learner, mode = "function")

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

        do.call(base_fun, call_args)
      }

    })

    assign(learner_name, new_fun, envir = parent.frame())
  }

  return(generated_learners)
}






#' Calculate Baseline Survival using Breslow Estimator
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (0=censored, 1=event).
#' @param risk_score Numeric vector of risk scores.
#' @param new.times Numeric vector of times at which to evaluate the baseline hazard.
#'
#' @return A numeric vector of cumulative baseline hazard #' computes \deqn{\hat{S}(t)}.
#' @keywords internal
#' @noRd
safe_breslow_step <- function(time, event, risk_score, new.times) {

  # 1. Sort data by time
  ord <- order(time, -event)
  time_sorted <- time[ord]
  event_sorted <- event[ord]

  # Safety: Clamp huge scores to avoid Inf.
  safe_score <- pmin(risk_score[ord], 700)
  risk_sorted <- exp(safe_score)

  # 2. Compute Risk Set Sums
  risk_pool <- rev(cumsum(rev(risk_sorted)))

  # 3. Compute Hazard contributions
  is_event <- event_sorted == 1
  if (!any(is_event)) return(rep(0, length(new.times)))

  event_t <- time_sorted[is_event]
  haz_contribution <- 1 / risk_pool[is_event]

  # 4. Aggregate ties
  unique_t <- unique(event_t)
  if (length(unique_t) < length(event_t)) {
    haz_contribution <- tapply(haz_contribution, event_t, sum)
    event_t <- as.numeric(names(haz_contribution))
  }

  # 5. Cumulative Hazard
  cum_haz <- cumsum(haz_contribution)

  # 6. Interpolate (Step Function)
  stats::approx(
    x = event_t,
    y = cum_haz,
    xout = new.times,
    method = "constant",
    yleft = 0,
    yright = max(cum_haz),
    f = 0,
    rule = 2,
    ties = mean
  )$y
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

