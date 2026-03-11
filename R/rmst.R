#' Calculate Restricted Mean Survival Time (RMST)
#' @param surv_matrix Matrix of survival probabilities (rows: patients, cols: time points)
#' @param times Vector of time points corresponding to the columns
#' @param tau The restriction time horizon
#' @return A vector of RMST values for each patient
#' @keywords internal
get_rmst <- function(surv_matrix, times, tau) {
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



#' Estimate Causal Restricted Mean Survival Time (RMST)
#'
#' Calculates the causal treatment effect based on the difference in Restricted Mean
#' Survival Time (RMST) between treatment and control groups up to a specific truncation time.
#'
#' @param fit A fitted \code{SuperSurv} ensemble object.
#' @param data A \code{data.frame} containing the patient covariates and the treatment assignment.
#' @param trt_col Character string. The exact name of the binary treatment indicator column in \code{data} (e.g., "treatment").
#' @param times Numeric vector of time points matching the prediction grid.
#' @param tau Numeric. A single truncation time limit up to which the RMST will be calculated.
#' @param verbose Logical; if \code{TRUE}, progress messages are shown.
#' @return A numeric value representing the estimated causal RMST difference (Treatment - Control).
#' @examples
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
#'   event.library = c("surv.coxph", "surv.glmnet"),
#'   cens.library = c("surv.coxph"),
#'   control = list(saveFitLibrary = TRUE)
#' )
#'
#' rmst_res <- estimate_marginal_rmst(
#'   fit = fit,
#'   data = dat,
#'   trt_col = "x4",
#'   times = new.times,
#'   tau = 100
#' )
#' rmst_res$ATE_RMST
#' @export
estimate_marginal_rmst <- function(fit, data, trt_col, times, tau,verbose = FALSE) {
  if(tau > max(times)) stop("tau cannot be greater than the maximum predicted time.")

  # Dynamically extract the exact variables the model was trained on
  model_vars <- fit$varNames

  # 1. Counterfactual: Everyone gets Exposure (A = 1)
  data_trt1 <- data
  data_trt1[[trt_col]] <- 1
  X_trt1 <- data_trt1[, model_vars, drop = FALSE]

  # Using the standardized 'newdata' argument
  pred_trt1_obj <- predict(fit, newdata = X_trt1, new.times = times)
  surv_trt1 <- pred_trt1_obj$event.predict
  rmst_1 <- get_rmst(surv_trt1, times, tau)

  # 2. Counterfactual: Everyone gets Control (A = 0)
  data_trt0 <- data
  data_trt0[[trt_col]] <- 0
  X_trt0 <- data_trt0[, model_vars, drop = FALSE]

  # Using the standardized 'newdata' argument
  pred_trt0_obj <- predict(fit, newdata = X_trt0, new.times = times)
  surv_trt0 <- pred_trt0_obj$event.predict
  rmst_0 <- get_rmst(surv_trt0, times, tau)

  # 3. Calculate the Difference
  ATE <- mean(rmst_1) - mean(rmst_0)

  res <- list(
    ATE_RMST = ATE,
    mean_RMST_Treated = mean(rmst_1),
    mean_RMST_Control = mean(rmst_0),
    tau = tau,
    patient_rmst_treated = rmst_1,
    patient_rmst_control = rmst_0
  )

  if (isTRUE(verbose)) {  message(sprintf( "Adjusted Delta RMST at tau=%s: %s time units", tau, round(ATE, 3))) }
  return(res)
}



#' Plot Causal RMST Difference Over Time
#'
#' Generates a curve showing how the causal Restricted Mean Survival Time (RMST) difference
#' between treatment groups evolves across a sequence of different truncation times.
#'
#' @param fit A fitted \code{SuperSurv} ensemble object.
#' @param data A \code{data.frame} containing the patient covariates and the treatment assignment.
#' @param trt_col Character string. The exact name of the binary treatment indicator column in \code{data}.
#' @param times Numeric vector of time points matching the prediction grid.
#' @param tau_seq Numeric vector. A sequence of truncation times (\code{tau}) to evaluate and plot.
#' @examples
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
#'   event.library = c("surv.coxph", "surv.glmnet"),
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
#' @return A \code{ggplot} object visualizing the causal RMST difference curve.
#' @export
plot_marginal_rmst_curve <- function(fit, data, trt_col, times, tau_seq) {
  requireNamespace("ggplot2", quietly = TRUE)
  results <- lapply(tau_seq, function(t) {
    res <- estimate_marginal_rmst(fit, data, trt_col, times, tau = t)
    data.frame(Tau = t, ATE = res$ATE_RMST)
  })
  res_df <- do.call(rbind, results)

  ggplot2::ggplot(res_df, ggplot2::aes(x = Tau, y = ATE)) +
    ggplot2::geom_line(color = "#e63946", linewidth = 1.2) +
    ggplot2::geom_point(color = "#1d3557", size = 3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Adjusted Effect over Time",
                  subtitle = "Difference in Restricted Mean Survival Time (Exposure 1 - 0)",
                  x = "Time Horizon (Tau)", y = expression(Delta ~ "RMST"))
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
#'   event.library = c("surv.coxph", "surv.glmnet"),
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
#'   tau = 350
#' )
#' @export
plot_rmst_vs_obs <- function(fit, data, time_col, event_col, times, tau) {
  requireNamespace("ggplot2", quietly = TRUE)
  model_vars <- fit$varNames
  X_obs <- data[, model_vars, drop = FALSE]

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
                  x = "Observed Survival Time", y = "Predicted RMST", color = "Status")
}
