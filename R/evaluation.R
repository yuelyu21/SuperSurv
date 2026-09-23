#' Marginal reverse Kaplan-Meier estimate for censoring
#' @noRd
.marginal_censoring_survival <- function(time, event) {
  grid <- sort(unique(time))
  survival <- cumprod(vapply(grid, function(value) {
    at_risk <- sum(time >= value)
    failures <- sum(time == value & event == 1)
    censorings <- sum(time == value & event == 0)
    denominator <- at_risk - failures
    if (censorings == 0) return(1)
    if (denominator <= 0) return(0)
    1 - censorings / denominator
  }, numeric(1L)))
  list(time = grid, survival = survival)
}


.validate_evaluation_inputs <- function(time, event, S_mat, times) {
  if (!is.numeric(time) || length(time) == 0L || any(!is.finite(time)) ||
      any(time < 0)) {
    stop("`time` must be a non-empty numeric vector of finite, non-negative values.",
         call. = FALSE)
  }
  if (!is.numeric(event) || length(event) != length(time) ||
      any(!is.finite(event)) || any(!event %in% c(0, 1))) {
    stop("`event` must be a numeric 0/1 vector with the same length as `time`.",
         call. = FALSE)
  }
  times <- .validate_time_grid(times, "times")
  if (!is.matrix(S_mat) && !is.data.frame(S_mat)) {
    stop("`S_mat` must be a numeric matrix with observations in rows and `times` in columns.",
         call. = FALSE)
  }
  S_mat <- as.matrix(S_mat)
  if (!is.numeric(S_mat)) {
    stop("`S_mat` must be a numeric matrix of survival probabilities.", call. = FALSE)
  }
  expected <- c(length(time), length(times))
  if (!identical(dim(S_mat), as.integer(expected))) {
    stop(
      "`S_mat` has dimensions ", paste(dim(S_mat), collapse = " x "),
      "; expected ", paste(expected, collapse = " x "),
      " to match `time` and `times`.",
      call. = FALSE
    )
  }
  if (any(!is.finite(S_mat)) || any(S_mat < 0 | S_mat > 1)) {
    stop("`S_mat` must contain finite survival probabilities in [0, 1].",
         call. = FALSE)
  }
  list(time = as.numeric(time), event = as.numeric(event),
       S_mat = S_mat, times = times)
}


.validate_integration_bounds <- function(tmin, tmax) {
  if (!is.numeric(tmin) || length(tmin) != 1L || !is.finite(tmin) ||
      !is.numeric(tmax) || length(tmax) != 1L || !is.finite(tmax) ||
      tmin > tmax) {
    stop("`tmin` and `tmax` must be finite numbers with `tmin <= tmax`.",
         call. = FALSE)
  }
  invisible(TRUE)
}


#' IPCW Brier Score and Integrated Brier Score (IBS)
#'
#' Calculates the Inverse Probability of Censoring Weighted (IPCW) Brier Score
#' over a grid of times, and computes the Integrated Brier Score (IBS) using
#' trapezoidal integration.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param S_mat A numeric matrix of predicted survival probabilities
#'   (rows = observations, columns = time points).
#' @param times Numeric vector of evaluation times matching the columns of \code{S_mat}.
#' @param tmin Numeric. Lower bound for IBS integration. Defaults to \code{min(times)}.
#' @param tmax Numeric. Upper bound for IBS integration. Defaults to \code{max(times)}.
#' @param ipcw_floor Positive numeric lower bound applied to the marginal
#'   Kaplan-Meier estimate of the censoring survival function before inversion.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{brier_scores}: A numeric vector of Brier scores at each time point.
#'   \item \code{ibs}: The Integrated Brier Score over the range `\code{tmin}, \code{tmax}`.
#'   \item \code{times}: The time grid used.
#' }
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:40, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:10, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.coxph(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' eval_brier(
#'   time = dat$duration[1:10],
#'   event = dat$event[1:10],
#'   S_mat = fit[["pred"]],
#'   times = times
#' )
#' @export
eval_brier <- function(time, event, S_mat, times, tmin = min(times),
                       tmax = max(times), ipcw_floor = 1e-6) {
  validated <- .validate_evaluation_inputs(time, event, S_mat, times)
  time <- validated$time
  event <- validated$event
  times <- validated$times
  S_mat <- validated$S_mat
  n <- length(time)
  T_len <- length(times)
  .validate_integration_bounds(tmin, tmax)
  if (length(ipcw_floor) != 1L || !is.finite(ipcw_floor) ||
      ipcw_floor <= 0 || ipcw_floor > 1) {
    stop("`ipcw_floor` must be one number in (0, 1].", call. = FALSE)
  }

  # Reverse Kaplan-Meier removes failures before censorings at tied times.
  fit_G <- .marginal_censoring_survival(time, event)
  t_fit <- fit_G$time
  s_fit <- fit_G$survival

  Gfun <- function(t, left_limit = FALSE) {
    pmax(.step_curve_at(t_fit, s_fit, t, left.limit = left_limit), ipcw_floor)
  }

  # Calculate Brier Score at each time t
  bs <- numeric(T_len)
  G_Ti_left <- Gfun(time, left_limit = TRUE)

  for (j in seq_len(T_len)) {
    t <- times[j]
    S_t <- S_mat[, j]

    Y <- as.integer(time > t)
    G_t <- Gfun(t)

    # Graf IPCW Weights
    w <- ifelse(time <= t & event == 1, 1 / G_Ti_left,
                ifelse(time > t, 1 / G_t, 0))

    bs[j] <- mean(w * (Y - S_t)^2)
  }

  # Calculate Integrated Brier Score (IBS) via Trapezoidal Rule
  valid_idx <- which(times >= tmin & times <= tmax)
  if (length(valid_idx) < 2) {
    ibs <- NA
  } else {
    t_sub <- times[valid_idx]
    b_sub <- bs[valid_idx]
    diffs <- diff(t_sub)
    heights <- (b_sub[-1] + b_sub[-length(b_sub)]) / 2
    ibs <- sum(diffs * heights) / (max(t_sub) - min(t_sub))
  }

  return(list(brier_scores = bs, ibs = ibs, times = times))
}


#' IPCW Log-Loss and Integrated IPCW Log-Loss
#'
#' Evaluates the two-term inverse-probability-of-censoring weighted log-loss.
#' By default, censoring survival is estimated using marginal reverse
#' Kaplan-Meier. Conditional or externally estimated censoring survival values
#' can instead be supplied explicitly. Probability clipping is applied only
#' inside logarithms.
#'
#' @inheritParams eval_brier
#' @param ipcw_floor Positive lower bound applied to censoring survival before
#'   inversion.
#' @param ipcw_cap Largest IPCW contribution. Use \code{Inf} for no cap.
#' @param eps Probability clipping constant in \code{(0, 0.5)}.
#' @param G_T_left Optional numeric vector containing subject-specific
#'   \eqn{G(T_i-\mid X_i)} values. Supply together with \code{G_times}.
#' @param G_times Optional numeric matrix with the dimensions of \code{S_mat}
#'   containing \eqn{G(t\mid X_i)}, or a vector aligned with \code{times} for a
#'   common censoring survival curve. Supply together with \code{G_T_left}.
#' @return A list containing \code{logloss_scores}, \code{integrated_logloss},
#'   \code{times}, and \code{diagnostics}. Diagnostics report the effective
#'   denominator threshold, intervention counts, weight quantiles, effective
#'   sample sizes, and the fraction of weighted loss contributed by capped
#'   rows, separately for failure and survivor terms.
#' @export
eval_logloss <- function(time, event, S_mat, times, tmin = min(times),
                         tmax = max(times), ipcw_floor = 1e-6,
                         eps = 1e-10, ipcw_cap = Inf,
                         G_T_left = NULL, G_times = NULL) {
  validated <- .validate_evaluation_inputs(time, event, S_mat, times)
  time <- validated$time
  event <- validated$event
  times <- validated$times
  S_mat <- validated$S_mat
  .validate_integration_bounds(tmin, tmax)
  if (length(ipcw_floor) != 1L || !is.finite(ipcw_floor) ||
      ipcw_floor <= 0 || ipcw_floor > 1) {
    stop("`ipcw_floor` must be one number in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(ipcw_cap) || length(ipcw_cap) != 1L ||
      is.na(ipcw_cap) || ipcw_cap < 1) {
    stop("`ipcw_cap` must be one number greater than or equal to 1, or Inf.",
         call. = FALSE)
  }
  if (length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("`eps` must be one number in (0, 0.5).", call. = FALSE)
  }

  supplied_G <- !is.null(G_T_left) || !is.null(G_times)
  if (xor(is.null(G_T_left), is.null(G_times))) {
    stop("Supply both `G_T_left` and `G_times`, or leave both NULL.",
         call. = FALSE)
  }
  if (!supplied_G) {
    fit_G <- .marginal_censoring_survival(time, event)
    Gfun <- function(x, left_limit = FALSE) {
      .step_curve_at(
        fit_G$time, fit_G$survival, x, left.limit = left_limit
      )
    }
    G_T_left <- Gfun(time, left_limit = TRUE)
    G_times <- matrix(
      rep(Gfun(times), each = length(time)),
      nrow = length(time), ncol = length(times)
    )
    censoring_source <- "marginal reverse Kaplan-Meier"
  } else {
    if (!is.numeric(G_T_left) || length(G_T_left) != length(time) ||
        any(!is.finite(G_T_left)) || any(G_T_left < 0 | G_T_left > 1)) {
      stop("`G_T_left` must contain one finite probability in [0, 1] per row.",
           call. = FALSE)
    }
    if (is.numeric(G_times) && is.null(dim(G_times)) &&
        length(G_times) == length(times)) {
      G_times <- matrix(
        rep(G_times, each = length(time)),
        nrow = length(time), ncol = length(times)
      )
    }
    if (!is.matrix(G_times) || !identical(dim(G_times), dim(S_mat)) ||
        !is.numeric(G_times) || any(!is.finite(G_times)) ||
        any(G_times < 0 | G_times > 1)) {
      stop(paste0(
        "`G_times` must be a finite probability matrix with the same ",
        "dimensions as `S_mat`, or a probability vector aligned with `times`."
      ), call. = FALSE)
    }
    censoring_source <- "user supplied"
  }

  n <- length(time)
  m <- length(times)
  failure_active <- outer(time, times, `<=`) & event == 1
  survivor_active <- outer(time, times, `>`)
  failure_denominator <- matrix(G_T_left, nrow = n, ncol = m)
  stabilized_failure <- .stabilize_ipcw(
    as.vector(failure_denominator), ipcw_floor, ipcw_cap
  )
  stabilized_survivor <- .stabilize_ipcw(
    as.vector(G_times), ipcw_floor, ipcw_cap
  )
  probability <- pmin(pmax(S_mat, eps), 1 - eps)
  failure_weight <- matrix(
    stabilized_failure$weight, nrow = n, ncol = m
  ) * failure_active
  survivor_weight <- matrix(
    stabilized_survivor$weight, nrow = n, ncol = m
  ) * survivor_active
  failure_contribution <- -failure_weight * log(1 - probability)
  survivor_contribution <- -survivor_weight * log(probability)
  contribution <- failure_contribution + survivor_contribution
  scores <- colMeans(contribution)

  valid <- which(times >= tmin & times <= tmax)
  integrated <- NA_real_
  if (length(valid) >= 2L) {
    selected_times <- times[valid]
    selected_scores <- scores[valid]
    integrated <- sum(
      diff(selected_times) *
        (selected_scores[-1L] + selected_scores[-length(selected_scores)]) / 2
    ) / diff(range(selected_times))
  }
  failure_active_vector <- as.vector(failure_active)
  survivor_active_vector <- as.vector(survivor_active)
  combined_active <- failure_active_vector | survivor_active_vector
  combined_stabilized <- list(
    probability = ifelse(
      failure_active_vector,
      stabilized_failure$probability,
      stabilized_survivor$probability
    ),
    raw.weight = ifelse(
      failure_active_vector,
      stabilized_failure$raw.weight,
      stabilized_survivor$raw.weight
    ),
    weight = ifelse(
      failure_active_vector,
      stabilized_failure$weight,
      stabilized_survivor$weight
    ),
    floored = ifelse(
      failure_active_vector,
      stabilized_failure$floored,
      stabilized_survivor$floored
    ),
    capped = ifelse(
      failure_active_vector,
      stabilized_failure$capped,
      stabilized_survivor$capped
    )
  )
  diagnostics <- list(
    censoring.source = censoring_source,
    ipcw.floor = ipcw_floor,
    ipcw.cap = ipcw_cap,
    effective.denominator.threshold = max(
      ipcw_floor, if (is.finite(ipcw_cap)) 1 / ipcw_cap else 0
    ),
    probability.clip.fraction = mean(S_mat < eps | S_mat > 1 - eps),
    failure.weights = .summarize_ipcw(
      stabilized_failure, failure_active_vector,
      as.vector(failure_contribution)
    ),
    survivor.weights = .summarize_ipcw(
      stabilized_survivor, survivor_active_vector,
      as.vector(survivor_contribution)
    ),
    combined.weights = .summarize_ipcw(
      combined_stabilized, combined_active, as.vector(contribution)
    )
  )

  list(
    logloss_scores = scores,
    integrated_logloss = integrated,
    times = times,
    diagnostics = diagnostics
  )
}






#' Calculate Concordance Index (Harrell's or Uno's)
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param S_mat A numeric matrix of predicted survival probabilities.
#' @param times Numeric vector of evaluation times matching the columns of \code{S_mat}.
#' @param eval_time Numeric. The specific time point at which to extract predictions.
#' @param method Character. Either "harrell" or "uno". Defaults to "uno".
#'
#' @return A numeric value representing the chosen C-index.
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:40, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#' X <- dat[, x_cols, drop = FALSE]
#' newX <- X[1:10, , drop = FALSE]
#' times <- seq(50, 150, by = 50)
#'
#' fit <- surv.coxph(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = newX,
#'   new.times = times,
#'   obsWeights = rep(1, nrow(dat)),
#'   id = NULL
#' )
#'
#' eval_cindex(
#'   time = dat$duration[1:10],
#'   event = dat$event[1:10],
#'   S_mat = fit[["pred"]],
#'   times = times,
#'   eval_time = 100,
#'   method = "uno"
#' )
#' @export
eval_cindex <- function(time, event, S_mat, times, eval_time, method = "uno") {
  validated <- .validate_evaluation_inputs(time, event, S_mat, times)
  time <- validated$time
  event <- validated$event
  S_mat <- validated$S_mat
  times <- validated$times
  method <- match.arg(method, c("uno", "harrell"))
  if (!is.numeric(eval_time) || length(eval_time) != 1L ||
      !is.finite(eval_time) || eval_time < min(times) || eval_time > max(times)) {
    stop("`eval_time` must be one finite number within the range of `times`.",
         call. = FALSE)
  }

  # Find the closest time column in the matrix
  t_idx <- which.min(abs(times - eval_time))

  # survival::concordance expects a "survival" score (higher = lives longer)
  # So we pass the survival probability DIRECTLY.
  surv_score <- S_mat[, t_idx]

  if (method == "uno") {
    # timewt = "n/G2" applies Uno's IPCW correction
    cfit <- survival::concordance(survival::Surv(time, event) ~ surv_score, timewt = "n/G2")
  } else {
    # Default Harrell's C
    cfit <- survival::concordance(survival::Surv(time, event) ~ surv_score)
  }

  return(cfit$concordance)
}






#' Time-Dependent AUC and Integrated AUC
#'
#' Evaluates the cumulative/dynamic time-dependent AUC and integrated AUC (iAUC)
#' using inverse probability of censoring weighting (IPCW).
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param S_mat A numeric matrix of predicted survival probabilities.
#' @param times Numeric vector of evaluation times matching the columns of \code{S_mat}.
#'
#' @return A list containing the \code{AUC_curve} at each time point, the
#'   \code{times}, and the integrated AUC \code{iAUC}.
#' @examples
#'  data("metabric", package = "SuperSurv")
#'  dat <- metabric[1:40, ]
#'  x_cols <- grep("^x", names(dat))[1:3]
#'  X <- dat[, x_cols, drop = FALSE]
#'  newX <- X[1:10, , drop = FALSE]
#'  times <- seq(50, 150, by = 50)
#'
#'  fit <- surv.coxph(
#'    time = dat$duration,
#'    event = dat$event,
#'    X = X,
#'    newdata = newX,
#'    new.times = times,
#'    obsWeights = rep(1, nrow(dat)),
#'    id = NULL
#'  )
#'
#'  eval_timeROC(
#'    time = dat$duration[1:10],
#'    event = dat$event[1:10],
#'    S_mat = fit[["pred"]],
#'    times = times
#'  )
#' @export
eval_timeROC <- function(time, event, S_mat, times) {
  validated <- .validate_evaluation_inputs(time, event, S_mat, times)
  time <- validated$time
  event <- validated$event
  S_mat <- validated$S_mat
  times <- validated$times
  n <- length(time)

  # Estimate censoring survival G(t) using reverse Kaplan-Meier.
  fit_G <- .marginal_censoring_survival(time, event)
  t_fit <- fit_G$time
  s_fit <- fit_G$survival
  eps_G <- 1e-6

  Gfun <- function(t, left_limit = FALSE) {
    pmax(.step_curve_at(t_fit, s_fit, t, left.limit = left_limit), eps_G)
  }

  auc_vals <- rep(NA_real_, length(times))

  for (j in seq_along(times)) {
    t0 <- times[j]

    # Use risk at the same evaluation time
    risk_score <- 1 - S_mat[, j]

    # Cumulative/dynamic definition
    is_case <- (time <= t0 & event == 1)
    is_ctrl <- (time > t0)

    if (sum(is_case) == 0 || sum(is_ctrl) == 0) {
      auc_vals[j] <- NA_real_
      next
    }

    # IPCW weights
    w_case <- rep(0, n)
    w_ctrl <- rep(0, n)

    w_case[is_case] <- 1 / Gfun(time[is_case], left_limit = TRUE)
    w_ctrl[is_ctrl] <- 1 / Gfun(t0)

    case_idx <- which(is_case)
    ctrl_idx <- which(is_ctrl)

    # Weighted pairwise AUC:
    # P(score_case > score_ctrl) + 0.5 P(tie)
    num <- 0
    den <- 0

    for (ii in case_idx) {
      for (jj in ctrl_idx) {
        wij <- w_case[ii] * w_ctrl[jj]
        den <- den + wij

        if (risk_score[ii] > risk_score[jj]) {
          num <- num + wij
        } else if (risk_score[ii] == risk_score[jj]) {
          num <- num + 0.5 * wij
        }
      }
    }

    auc_vals[j] <- if (den > 0) num / den else NA_real_
  }

  valid_idx <- which(!is.na(auc_vals))
  if (length(valid_idx) < 2) {
    iauc <- NA_real_
  } else {
    t_sub <- times[valid_idx]
    a_sub <- auc_vals[valid_idx]
    diffs <- diff(t_sub)
    heights <- (a_sub[-1] + a_sub[-length(a_sub)]) / 2
    iauc <- sum(diffs * heights) / (max(t_sub) - min(t_sub))
  }

  list(
    AUC_curve = auc_vals,
    times = times,
    iAUC = iauc
  )
}





#' Evaluate survival predictions across models and times
#'
#' Computes a common set of numerical benchmark results for a fitted
#' \code{SuperSurv} object or supported standalone learner. The censoring
#' distribution used by the Brier score and time-dependent AUC is estimated
#' marginally by Kaplan-Meier. For conditional censoring models, resampling,
#' inference, or formal model comparisons, use a specialist evaluator such as
#' \code{riskRegression::Score()}.
#'
#' @inheritParams eval_summary
#' @return A list of class \code{"SuperSurv_benchmark"} containing
#'   \code{summary}, a model-level table; \code{by_time}, a time-specific table;
#'   and the prediction grid and risk horizon.
#' @export
eval_benchmark <- function(object, newdata, time, event, eval_times,
                           risk_time = stats::median(eval_times),
                           verbose = FALSE) {
  if (!inherits(object, "SuperSurv") && !is.object(object)) {
    stop("`object` must be a fitted 'SuperSurv' object or a supported fitted survival learner.",
         call. = FALSE)
  }
  if (inherits(object, "SuperSurv")) {
    .validate_SuperSurv_object(object)
    newdata <- .validate_prediction_newdata(newdata, object)
  } else {
    if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
      stop("`newdata` must be a non-empty data frame.", call. = FALSE)
    }
    .validate_data_frame_columns(newdata, "newdata")
  }
  eval_times <- .validate_time_grid(eval_times, "eval_times")
  if (!is.numeric(time) || length(time) != nrow(newdata) ||
      any(!is.finite(time)) || any(time < 0)) {
    stop("`time` must contain one finite, non-negative value per row of `newdata`.",
         call. = FALSE)
  }
  if (!is.numeric(event) || length(event) != nrow(newdata) ||
      any(!is.finite(event)) || any(!event %in% c(0, 1))) {
    stop("`event` must contain one numeric 0/1 value per row of `newdata`.",
         call. = FALSE)
  }
  if (!is.numeric(risk_time) || length(risk_time) != 1L ||
      !is.finite(risk_time) || risk_time < min(eval_times) ||
      risk_time > max(eval_times)) {
    stop("`risk_time` must be one finite number within the range of `eval_times`.",
         call. = FALSE)
  }
  .validate_scalar_logical(verbose, "verbose")
  if (isTRUE(verbose)) message("Generating predictions on test data...")
  preds <- tryCatch(
    predict(object, newdata = newdata, new.times = eval_times),
    error = function(error) {
      stop("Could not generate predictions from `object`: ",
           conditionMessage(error), call. = FALSE)
    }
  )

  prediction_matrices <- list()
  if (is.list(preds) && !is.null(preds$event.predict)) {
    prediction_matrices[["SuperSurv_Ensemble"]] <- as.matrix(preds$event.predict)
    library_array <- preds$event.library.predict
    k_models <- dim(library_array)[3]
    model_names <- dimnames(library_array)[[3]]
    if (is.null(model_names) || length(model_names) != k_models) {
      model_names <- object$event.libraryNames
    }
    if (is.data.frame(model_names)) model_names <- model_names[[1L]]
    if (is.null(model_names) || length(model_names) != k_models) {
      model_names <- paste0("Base_Learner_", seq_len(k_models))
    }
    model_names <- make.unique(as.character(model_names))
    for (i in seq_len(k_models)) {
      prediction_matrices[[model_names[i]]] <- library_array[, , i, drop = TRUE]
    }
  } else if (is.matrix(preds)) {
    model_name <- class(object)[1L]
    if (is.null(model_name) || identical(model_name, "matrix")) {
      model_name <- "Standalone_Model"
    }
    prediction_matrices[[model_name]] <- preds
  } else {
    stop("Unrecognized prediction format from `object`.", call. = FALSE)
  }

  evaluate_one <- function(model_name, S_mat) {
    brier <- eval_brier(time, event, S_mat, eval_times)
    logloss <- eval_logloss(time, event, S_mat, eval_times)
    auc <- suppressMessages(eval_timeROC(time, event, S_mat, eval_times))
    cindex_curve <- vapply(eval_times, function(t) {
      eval_cindex(time, event, S_mat, eval_times, eval_time = t, method = "uno")
    }, numeric(1L))
    uno_c <- eval_cindex(
      time, event, S_mat, eval_times, eval_time = risk_time, method = "uno"
    )
    list(
      summary = data.frame(
        Model = model_name, IBS = brier$ibs,
        IPCW_LogLoss = logloss$integrated_logloss,
        Uno_C = uno_c, iAUC = auc$iAUC,
        stringsAsFactors = FALSE
      ),
      by_time = data.frame(
        Time = eval_times, Model = model_name, Brier = brier$brier_scores,
        LogLoss = logloss$logloss_scores,
        CD_AUC = auc$AUC_curve, C_Index = cindex_curve,
        stringsAsFactors = FALSE
      )
    )
  }

  evaluated <- Map(evaluate_one, names(prediction_matrices), prediction_matrices)
  result <- list(
    summary = do.call(rbind, lapply(evaluated, `[[`, "summary")),
    by_time = do.call(rbind, lapply(evaluated, `[[`, "by_time")),
    eval_times = as.numeric(eval_times),
    risk_time = risk_time
  )
  rownames(result$summary) <- NULL
  rownames(result$by_time) <- NULL
  class(result) <- "SuperSurv_benchmark"
  result
}


#' Evaluate SuperSurv predictions on test data
#'
#' Computes the integrated Brier score (IBS), Uno C-index, and integrated area
#' under the curve (iAUC) for the SuperSurv ensemble and all individual base
#' learners.
#'
#' @param object A fitted \code{SuperSurv} object.
#' @param newdata A data.frame of test covariates.
#' @param time Numeric vector of observed follow-up times for the test set.
#' @param event Numeric vector of event indicators for the test set.
#' @param eval_times Numeric vector of times at which to evaluate survival
#'   predictions.
#' @param risk_time Numeric. The specific time horizon used when extracting risk
#'   scores for Uno C-index. Defaults to the median of \code{eval_times}.
#' @param verbose Logical; if \code{TRUE}, progress messages are shown.
#'
#' @return A data frame of integrated benchmark metrics for the ensemble and
#'   base learners. Use \code{eval_benchmark()} to also obtain time-specific
#'   numerical results.
#'
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:80, ]
#' x_cols <- grep("^x", names(dat))[1:3]
#'
#' fit <- SuperSurv(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = dat[, x_cols, drop = FALSE],
#'   new.times = seq(50, 200, by = 50),
#'   event.library = c("surv.coxph", "surv.km"),
#'   cens.library = c("surv.coxph", "surv.km")
#' )
#'
#' res <- eval_summary(
#'   object = fit,
#'   newdata = dat[, x_cols, drop = FALSE],
#'   time = dat$duration,
#'   event = dat$event,
#'   eval_times = seq(50, 200, by = 50)
#' )
#'
#' res
#'
#' @export
eval_summary <- function(object, newdata, time, event, eval_times,
                         risk_time = stats::median(eval_times),
                         verbose = FALSE) {
  benchmark <- eval_benchmark(
    object = object, newdata = newdata, time = time, event = event,
    eval_times = eval_times, risk_time = risk_time, verbose = verbose
  )
  results <- benchmark$summary[c("Model", "IBS", "Uno_C", "iAUC")]
  results[c("IBS", "Uno_C", "iAUC")] <- lapply(
    results[c("IBS", "Uno_C", "iAUC")], round, digits = 4
  )

  attr(results, "eval_times") <- eval_times
  attr(results, "risk_time") <- risk_time
  attr(results, "by_time") <- benchmark$by_time
  class(results) <- c("SuperSurv_eval", class(results))

  results
}

