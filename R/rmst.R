#' Calculate Restricted Mean Survival Time (RMST)
#' @param surv_matrix Matrix of survival probabilities (rows: patients, cols: time points)
#' @param times Vector of time points corresponding to the columns
#' @param tau The restriction time horizon
#' @return A vector of RMST values for each patient
#' @keywords internal
get_rmst <- function(surv_matrix, times, tau) {
  times <- .validate_time_grid(times, "times")
  if (!is.matrix(surv_matrix) && !is.data.frame(surv_matrix)) {
    stop("`surv_matrix` must be a numeric matrix with observations in rows and `times` in columns.",
         call. = FALSE)
  }
  surv_matrix <- as.matrix(surv_matrix)
  if (!is.numeric(surv_matrix) || ncol(surv_matrix) != length(times) ||
      nrow(surv_matrix) == 0L) {
    stop("`surv_matrix` must be a non-empty numeric matrix with one column per value in `times`.",
         call. = FALSE)
  }
  if (any(!is.finite(surv_matrix)) || any(surv_matrix < 0 | surv_matrix > 1)) {
    stop("`surv_matrix` must contain finite survival probabilities in [0, 1].",
         call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau < 0 || tau > max(times)) {
    stop("`tau` must be one finite, non-negative number no greater than `max(times)`.",
         call. = FALSE)
  }
  surv_matrix <- as.matrix(surv_matrix)
  ord <- order(times)
  times <- times[ord]
  surv_matrix <- surv_matrix[, ord, drop = FALSE]

  if (tau <= 0) return(rep(0, nrow(surv_matrix)))
  if (tau <= min(times)) return(rep(tau, nrow(surv_matrix)))  # assume S(t)=1 before first grid

  valid_idx <- which(times <= tau)
  t_valid <- times[valid_idx]
  S_valid <- surv_matrix[, valid_idx, drop = FALSE]

  t_calc <- c(0, t_valid, tau)
  dt <- diff(t_calc)

  S_calc <- cbind(1, S_valid, S_valid[, ncol(S_valid), drop = FALSE])
  S_left <- as.matrix(S_calc[, -ncol(S_calc), drop = FALSE])

  rmst <- as.vector(S_left %*% matrix(dt, ncol = 1))
  rmst
}




#' Estimate an Adjusted Marginal RMST Contrast
#'
#' Computes a covariate-adjusted marginal contrast on the restricted mean
#' survival time (RMST) scale using standardization (g-computation) based on a
#' fitted \code{SuperSurv} model.
#'
#' For a binary grouping variable \code{trt_col}, the function predicts
#' counterfactual survival curves under \code{A = 1} and \code{A = 0} for every
#' individual in the supplied dataset, integrates each curve up to the
#' restriction time \code{tau}, and averages the resulting individual-level
#' RMST differences. The resulting contrast is generally interpreted as an
#' adjusted marginal contrast. When \code{trt_col} corresponds to a manipulable
#' intervention and additional identification assumptions hold, the same
#' standardized procedure may also support a causal interpretation.
#'
#' @param fit A fitted object of class \code{"SuperSurv"}.
#' @param data A \code{data.frame} containing the covariates used for
#'   standardization, including the binary grouping variable specified by
#'   \code{trt_col}.
#' @param trt_col Character string giving the name of the binary grouping
#'   variable in \code{data}. The variable is set to 1 and 0, respectively, to
#'   generate the two standardized prediction regimes.
#' @param times Numeric vector of prediction time points corresponding to the
#'   evaluation grid used for survival prediction.
#' @param tau Numeric scalar giving the restriction horizon for RMST. Must not
#'   exceed \code{max(times)}.
#' @param inference Deprecated compatibility argument. Formal inference is not
#'   provided; it requires a validated procedure accounting for model-estimation
#'   uncertainty. Must remain \code{FALSE}.
#' @param B,seed,ci_level Deprecated compatibility arguments; ignored when
#'   \code{inference = FALSE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{ATE_RMST}: The estimated adjusted marginal RMST contrast
#'   \eqn{\widehat{\Delta}_{RMST}(\tau)}.
#'   \item \code{mean_RMST_Treated}: The average predicted RMST under
#'   \code{A = 1}.
#'   \item \code{mean_RMST_Control}: The average predicted RMST under
#'   \code{A = 0}.
#'   \item \code{tau}: The restriction horizon used for integration.
#'   \item \code{patient_rmst_treated}: Vector of individual-level predicted
#'   RMST values under \code{A = 1}.
#'   \item \code{patient_rmst_control}: Vector of individual-level predicted
#'   RMST values under \code{A = 0}.
#'   \item \code{patient_delta_rmst}: Vector of individual-level predicted RMST
#'   contrasts.
#'   \item \code{inference}: Always \code{FALSE}; retained for compatibility.
#' }
#'
#' @details
#' The function uses the empirical distribution of the observed covariates in
#' \code{data} as the standardization distribution. RMST is evaluated
#' numerically from the predicted survival matrix using a left Riemann sum over
#' the supplied grid \code{times}.
#'
#' The returned contrast is a model-based point estimate. Formal uncertainty
#' quantification would need to account for nuisance estimation, tuning,
#' cross-validation, learner fitting, and ensemble-weight estimation, for
#' example through a validated full-pipeline refitting procedure.
#'
#' @examples
#' \dontrun{
#' data("metabric", package = "SuperSurv")
#' x_cols <- grep("^x", names(metabric), value = TRUE)
#' X <- metabric[, x_cols]
#' new.times <- seq(10, 150, by = 10)
#'
#' fit <- SuperSurv(
#'   time = metabric$duration,
#'   event = metabric$event,
#'   X = X,
#'   newdata = X,
#'   new.times = new.times,
#'   event.library = c("surv.coxph", "surv.rfsrc"),
#'   cens.library = c("surv.coxph"),
#'   control = list(saveFitLibrary = TRUE),
#'   nFolds = 3
#' )
#'
#' rmst_res <- estimate_marginal_rmst(
#'   fit = fit,
#'   data = metabric,
#'   trt_col = "x4",
#'   times = new.times,
#'   tau = 100
#' )
#'
#' rmst_res$ATE_RMST
#' }
#'
#' @export
estimate_marginal_rmst <- function(fit, data, trt_col, times, tau,
                                 inference = FALSE, B = 200,
                                 seed = NULL, ci_level = 0.95) {
  .validate_SuperSurv_object(fit, "fit")

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(trt_col) || length(trt_col) != 1L || is.na(trt_col) ||
      !nzchar(trt_col)) {
    stop("`trt_col` must be one non-empty column name.", call. = FALSE)
  }
  .validate_scalar_logical(inference, "inference")
  if (isTRUE(inference)) {
    stop(
      "Perturbation-based RMST inference has been removed because it does not account for model fitting, tuning, cross-validation, or ensemble selection. Use the returned RMST point contrast; this package does not provide a validated inference procedure accounting for these sources of uncertainty.",
      call. = FALSE
    )
  }
  if (!trt_col %in% names(data)) {
    stop("Treatment column '", trt_col, "' was not found in `data`.", call. = FALSE)
  }
  treatment <- data[[trt_col]]
  if (anyNA(treatment) || !all(treatment %in% c(0, 1))) {
    stop("`data[[trt_col]]` must contain only 0 and 1 with no missing values.",
         call. = FALSE)
  }
  times <- .validate_time_grid(times, "times")
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau < 0 || tau > max(times)) {
    stop("`tau` must be one finite, non-negative number no greater than `max(times)`.",
         call. = FALSE)
  }
  # Use the exact variables the model was trained on.
  model_vars <- training_variables(fit)

  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following training variables are missing from 'data': ",
         paste(missing_vars, collapse = ", "))
  }
  .validate_data_frame_columns(data[, model_vars, drop = FALSE], "data")

  # -----------------------------
  # Counterfactual predictions: A = 1
  # -----------------------------
  data_trt1 <- data
  data_trt1[[trt_col]] <- 1
  X_trt1 <- data_trt1[, model_vars, drop = FALSE]

  pred_trt1_obj <- predict(fit, newdata = X_trt1, new.times = times)
  surv_trt1 <- pred_trt1_obj$event.predict

  rmst_1 <- get_rmst(surv_matrix = surv_trt1, times = times, tau = tau)

  # -----------------------------
  # Counterfactual predictions: A = 0
  # -----------------------------
  data_trt0 <- data
  data_trt0[[trt_col]] <- 0
  X_trt0 <- data_trt0[, model_vars, drop = FALSE]

  pred_trt0_obj <- predict(fit, newdata = X_trt0, new.times = times)
  surv_trt0 <- pred_trt0_obj$event.predict

  rmst_0 <- get_rmst(surv_matrix = surv_trt0, times = times, tau = tau)

  # -----------------------------
  # Point estimate
  # -----------------------------
  delta_i <- rmst_1 - rmst_0
  ATE <- mean(delta_i)

  res <- list(
    ATE_RMST = ATE,
    mean_RMST_Treated = mean(rmst_1),
    mean_RMST_Control = mean(rmst_0),
    tau = tau,
    patient_rmst_treated = rmst_1,
    patient_rmst_control = rmst_0,
    patient_delta_rmst = delta_i,
    inference = FALSE
  )

  msg <- sprintf("Adjusted Delta RMST at tau = %s: %s time units",
                 tau, round(ATE, 3))

  message(msg)

  return(res)
}








#' Plot Adjusted Marginal RMST Contrast Over Time
#'
#' Generates a curve showing how the adjusted marginal restricted mean survival
#' time (RMST) contrast evolves across a sequence of restriction times.
#'
#' @param fit A fitted \code{SuperSurv} ensemble object.
#' @param data A \code{data.frame} containing the covariates and the binary grouping variable.
#' @param trt_col Character string. The exact name of the binary grouping variable in \code{data}.
#' @param times Numeric vector of time points matching the prediction grid.
#' @param tau_seq Numeric vector. A sequence of restriction times (\code{tau}) to evaluate and plot.
#' @param inference Deprecated compatibility argument. Must remain \code{FALSE}.
#' @param B,seed,ci_level Deprecated compatibility arguments; ignored.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:80, ]
#' x_cols <- grep("^x", names(dat), value = TRUE)[1:5]
#' X <- dat[, x_cols, drop = FALSE]
#' new.times <- seq(20, 120, by = 20)
#'
#' fit <- SuperSurv(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = X,
#'   new.times = new.times,
#'   event.library = c("surv.coxph"),
#'   cens.library = c("surv.coxph"),
#'   control = list(saveFitLibrary = TRUE)
#' )
#'
#' tau_grid <- seq(40, 120, by = 20)
#' plot_marginal_rmst_curve(
#'   fit = fit,
#'   data = dat,
#'   trt_col = "x4",
#'   times = new.times,
#'   tau_seq = tau_grid
#' )
#' }
#'
#' @return A \code{ggplot} object visualizing the adjusted marginal RMST contrast curve.
#' @export
plot_marginal_rmst_curve <- function(fit, data, trt_col, times, tau_seq,
                                     inference = FALSE, B = 200,
                                     seed = NULL, ci_level = 0.95) {
  .require_optional_packages("ggplot2", "plot_marginal_rmst_curve()")
  .validate_scalar_logical(inference, "inference")
  if (isTRUE(inference)) {
    stop(
      "Perturbation-based RMST confidence intervals have been removed because they do not account for the fitted learning and ensemble-selection pipeline.",
      call. = FALSE
    )
  }
  times <- .validate_time_grid(times, "times")
  tau_seq <- .validate_time_grid(tau_seq, "tau_seq")
  if (max(tau_seq) > max(times)) {
    stop("Every value in `tau_seq` must be no greater than `max(times)`.",
         call. = FALSE)
  }

  results <- lapply(seq_along(tau_seq), function(i) {
    t <- tau_seq[i]

    res <- estimate_marginal_rmst(
      fit = fit,
      data = data,
      trt_col = trt_col,
      times = times,
      tau = t
    )

    out <- data.frame(
      Tau = t,
      Delta_RMST = if (!is.null(res$Delta_RMST)) res$Delta_RMST else res$ATE_RMST,
      stringsAsFactors = FALSE
    )

    out
  })

  res_df <- do.call(rbind, results)

  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = Tau, y = Delta_RMST))

  p <- p +
    ggplot2::geom_line(color = "#e63946", linewidth = 1) +
    ggplot2::geom_point(color = "#1d3557", size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Adjusted Marginal RMST Contrast Over Time",
      subtitle = "Model-based difference in restricted mean survival time",
      x = "Restriction Time (Tau)",
      y = expression(Delta ~ "RMST")
    ) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", size = 16),
      plot.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.title = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 16)
    )

  return(p)
}


#' Plot Predicted RMST vs. Observed Survival Times
#'
#' Evaluates the calibration of the causal RMST estimator by plotting the model's
#' predicted RMST for each patient against their actual observed follow-up time.
#'
#' @param fit A fitted \code{SuperSurv} ensemble object.
#' @param data A \code{data.frame} containing the patient covariates, times, and events.
#' @param time_col Character string. The exact name of the observed follow-up time column in \code{data}.
#' @param event_col Character string. The exact name of the event indicator column in \code{data} (e.g., 1 for event, 0 for censored).
#' @param times Numeric vector of time points matching the prediction grid.
#' @param tau Numeric. A single truncation time limit up to which the RMST is calculated.
#'
#' @return A \code{ggplot} object comparing predicted RMST to observed outcomes.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:80, ]
#' x_cols <- grep("^x", names(dat))[1:5]
#' X <- dat[, x_cols, drop = FALSE]
#' new.times <- seq(20, 120, by = 20)
#'
#' fit <- SuperSurv(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   newdata = X,
#'   new.times = new.times,
#'   event.library = c("surv.coxph"),
#'   cens.library = c("surv.coxph"),
#'   control = list(saveFitLibrary = TRUE)
#' )
#'
#' plot_rmst_vs_obs(
#'   fit = fit,
#'   data = dat,
#'   time_col = "duration",
#'   event_col = "event",
#'   times = new.times,
#'   tau = 100
#' )
#' }
#' @export
plot_rmst_vs_obs <- function(fit, data, time_col, event_col, times, tau) {
  .require_optional_packages("ggplot2", "plot_rmst_vs_obs()")
  .validate_SuperSurv_object(fit, "fit")
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  }
  for (argument in c("time_col", "event_col")) {
    value <- get(argument)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("`", argument, "` must be one non-empty column name.", call. = FALSE)
    }
    if (!value %in% names(data)) {
      stop("Column '", value, "' supplied through `", argument,
           "` was not found in `data`.", call. = FALSE)
    }
  }
  observed_time <- data[[time_col]]
  observed_event <- data[[event_col]]
  if (!is.numeric(observed_time) || any(!is.finite(observed_time)) ||
      any(observed_time < 0)) {
    stop("`data[[time_col]]` must contain finite, non-negative numeric values.",
         call. = FALSE)
  }
  if (!is.numeric(observed_event) || any(!is.finite(observed_event)) ||
      any(!observed_event %in% c(0, 1))) {
    stop("`data[[event_col]]` must contain only numeric 0/1 values.",
         call. = FALSE)
  }
  times <- .validate_time_grid(times, "times")
  model_vars <- training_variables(fit)
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars)) {
    stop("`data` is missing training feature(s): ",
         paste(missing_vars, collapse = ", "), ".", call. = FALSE)
  }
  X_obs <- data[, model_vars, drop = FALSE]
  .validate_data_frame_columns(X_obs, "data")

  # Using the standardized 'newdata' argument
  pred_obs_obj <- predict(fit, newdata = X_obs, new.times = times)
  surv_obs <- pred_obs_obj$event.predict
  rmst_obs <- get_rmst(surv_obs, times, tau)

  plot_df <- data.frame(
    Observed_Time = data[[time_col]],
    Predicted_RMST = rmst_obs,
    Event = as.factor(data[[event_col]])
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Observed_Time, y = Predicted_RMST, color = Event)) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = tau, linetype = "dotted", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("0" = "#457b9d", "1" = "#e63946"),
                                labels = c("0" = "Censored", "1" = "Event")) +
    ggplot2::labs(title = sprintf("Predicted RMST vs. Observed Time (Tau = %s)", tau),
                  x = "Observed Survival Time", y = "Predicted RMST", color = "Status") +
   ggplot2::theme(
    legend.position = "right",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "bold",size = 16),
    plot.title = ggplot2::element_text(size = 20, face = "bold"),
    axis.title = ggplot2::element_text(size = 18),
    axis.text.x = ggplot2::element_text(size = 16)
  )
}
