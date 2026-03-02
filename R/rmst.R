#' Calculate Restricted Mean Survival Time (RMST)
#' @param surv_matrix Matrix of survival probabilities (rows: patients, cols: time points)
#' @param times Vector of time points corresponding to the columns
#' @param tau The restriction time horizon
#' @return A vector of RMST values for each patient
#' @keywords internal
get_rmst <- function(surv_matrix, times, tau) {

  # 1. Force the input into a strict numeric matrix
  surv_matrix <- as.matrix(surv_matrix)

  valid_idx <- which(times <= tau)
  t_valid <- times[valid_idx]

  # 2. Extract the valid columns, keeping it as a matrix
  S_valid <- surv_matrix[, valid_idx, drop = FALSE]

  t_calc <- c(0, t_valid, tau)

  # 3. Calculate time intervals and force it to be a column matrix
  dt <- diff(t_calc)
  dt_matrix <- matrix(dt, ncol = 1)

  # 4. Construct the step function probability matrix
  # At time 0, survival is exactly 1
  S_calc <- cbind(1, S_valid, S_valid[, ncol(S_valid), drop = FALSE])

  # Remove the last column to match the number of time intervals
  S_final <- as.matrix(S_calc[, -ncol(S_calc), drop = FALSE])

  # 5. Matrix multiplication: (N x time_steps) %*% (time_steps x 1)
  rmst <- as.vector(S_final %*% dt_matrix)

  return(rmst)
}





#' Estimate Causal Treatment Effect via RMST Difference
#' @export
estimate_causal_rmst <- function(fit, data, trt_col, times, tau) {

  if(tau > max(times)) stop("tau cannot be greater than the maximum predicted time.")

  # 1. Counterfactual: Everyone gets Treatment (A = 1)
  data_trt1 <- data
  data_trt1[[trt_col]] <- 1
  X_trt1 <- data_trt1[, grep("^x|A", names(data_trt1))]

  # Predict using 'newdata' and 'new.times'
  pred_trt1_obj <- predict(fit, newdata = X_trt1, new.times = times)
  surv_trt1 <- pred_trt1_obj$event.SL.predict

  rmst_1 <- get_rmst(surv_trt1, times, tau)

  # 2. Counterfactual: Everyone gets Control (A = 0)
  data_trt0 <- data
  data_trt0[[trt_col]] <- 0
  X_trt0 <- data_trt0[, grep("^x|A", names(data_trt0))]

  # Predict using 'newdata' and 'new.times'
  pred_trt0_obj <- predict(fit, newdata = X_trt0, new.times = times)
  surv_trt0 <- pred_trt0_obj$event.SL.predict

  rmst_0 <- get_rmst(surv_trt0, times, tau)

  # 3. Calculate Causal Difference (Average Treatment Effect)
  ATE <- mean(rmst_1) - mean(rmst_0)

  res <- list(
    ATE_RMST = ATE,
    mean_RMST_Treated = mean(rmst_1),
    mean_RMST_Control = mean(rmst_0),
    tau = tau,
    patient_rmst_treated = rmst_1,
    patient_rmst_control = rmst_0
  )

  message(sprintf("Causal ATE (Delta RMST at tau=%s): %s time units", tau, round(ATE, 3)))
  return(res)
}




#' Plot Causal Treatment Effect over Time (Delta RMST)
#' @export
plot_causal_rmst_curve <- function(fit, data, trt_col, times, tau_seq) {
  requireNamespace("ggplot2", quietly = TRUE)

  # Calculate ATE for every tau in the sequence
  results <- lapply(tau_seq, function(t) {
    res <- estimate_causal_rmst(fit, data, trt_col, times, tau = t)
    data.frame(Tau = t, ATE = res$ATE_RMST,
               Treated = res$mean_RMST_Treated,
               Control = res$mean_RMST_Control)
  })

  res_df <- do.call(rbind, results)

  # Plot the Delta RMST
  ggplot2::ggplot(res_df, ggplot2::aes(x = Tau, y = ATE)) +
    ggplot2::geom_line(color = "#e63946", size = 1.2) +
    ggplot2::geom_point(color = "#1d3557", size = 3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Causal Treatment Effect over Time",
      subtitle = "Difference in Restricted Mean Survival Time (Treated - Control)",
      x = "Time Horizon (Tau)",
      y = expression(Delta ~ "RMST")
    )
}



#' Plot Predicted RMST vs Observed Survival Time
#' @export
plot_rmst_vs_obs <- function(fit, data, time_col, event_col, times, tau) {
  requireNamespace("ggplot2", quietly = TRUE)

  X_obs <- data[, grep("^x|A", names(data))]

  # Predict using 'newdata' and 'new.times'
  pred_obs_obj <- predict(fit, newdata = X_obs, new.times = times)
  surv_obs <- pred_obs_obj$event.SL.predict

  # Calculate RMST
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
    ggplot2::labs(
      title = sprintf("Predicted RMST vs. Observed Time (Tau = %s)", tau),
      x = "Observed Survival Time",
      y = "Predicted RMST",
      color = "Status"
    )
}
