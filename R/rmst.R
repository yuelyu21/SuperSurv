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
#' If \code{inference = TRUE}, the function additionally performs a
#' perturbation-based inference procedure conditional on the fitted
#' \code{SuperSurv} model. In this implementation, the fitted learner library,
#' hyperparameters, base learners, and ensemble weights are held fixed, and
#' random positive weights are applied to the individual-level RMST contrasts to
#' estimate a perturbation-based standard error, Wald-type confidence interval,
#' and p-value.
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
#' @param inference Logical; if \code{TRUE}, compute perturbation-based standard
#'   errors, confidence intervals, and a Wald-type p-value. Defaults to
#'   \code{FALSE}.
#' @param B Integer giving the number of perturbation replicates when
#'   \code{inference = TRUE}. Defaults to \code{200}.
#' @param seed Optional integer seed for reproducibility of the perturbation
#'   procedure.
#' @param ci_level Numeric scalar in \code{(0, 1)} specifying the confidence
#'   level for the Wald-type confidence interval. Defaults to \code{0.95}.
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
#'   \item \code{inference}: Logical indicator for whether perturbation-based
#'   inference was requested.
#'   \item \code{B}: Number of perturbation replicates used when
#'   \code{inference = TRUE}; otherwise \code{NULL}.
#'   \item \code{SE_RMST}: Perturbation-based standard error of the RMST
#'   contrast; otherwise \code{NULL}.
#'   \item \code{CI_RMST}: Wald-type confidence interval for the RMST contrast;
#'   otherwise \code{NULL}.
#'   \item \code{z_value}: Wald-type test statistic; otherwise \code{NULL}.
#'   \item \code{p_value}: Two-sided Wald-type p-value; otherwise \code{NULL}.
#'   \item \code{perturb_reps}: Vector of perturbation replicate estimates;
#'   otherwise \code{NULL}.
#' }
#'
#' @details
#' The function uses the empirical distribution of the observed covariates in
#' \code{data} as the standardization distribution. RMST is evaluated
#' numerically from the predicted survival matrix using a left Riemann sum over
#' the supplied grid \code{times}.
#'
#' The perturbation-based inference implemented here is conditional on the
#' fitted \code{SuperSurv} model. It does not re-tune hyperparameters, reselect
#' the learner library, or refit the base learners under each perturbation.
#' Instead, it perturbs the aggregation of the individual-level predicted RMST
#' contrasts. This yields a lightweight uncertainty quantification procedure for
#' the standardized RMST contrast given the final fitted ensemble.
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
#'   tau = 100,
#'   inference = TRUE,
#'   B = 200,
#'   seed = 123
#' )
#'
#' rmst_res$ATE_RMST
#' rmst_res$SE_RMST
#' rmst_res$CI_RMST
#' format.pval(rmst_res$p_value, digits = 3, eps = 1e-16)
#' }
#'
#' @export
estimate_marginal_rmst <- function(fit, data, trt_col, times, tau,
                                 inference = FALSE, B = 200,
                                 seed = NULL, ci_level = 0.95) {

  if (!inherits(fit, "SuperSurv")) {
    stop("'fit' must be a fitted 'SuperSurv' object.")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }

  if (!trt_col %in% names(data)) {
    stop("'", trt_col, "' not found in 'data'.")
  }

  if (tau > max(times)) {
    stop("'tau' cannot be greater than the maximum value in 'times'.")
  }

  if (!is.logical(inference) || length(inference) != 1) {
    stop("'inference' must be TRUE or FALSE.")
  }

  if (!is.numeric(B) || length(B) != 1 || B <= 1) {
    stop("'B' must be a single integer greater than 1.")
  }

  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a single number in (0, 1).")
  }

  # Use the exact variables the model was trained on
  model_vars <- fit$varNames
  if (is.null(model_vars)) {
    stop("Could not find 'fit$varNames'. The fitted SuperSurv object must store training variable names.")
  }

  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following training variables are missing from 'data': ",
         paste(missing_vars, collapse = ", "))
  }

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

  # -----------------------------
  # Optional perturbation inference
  # -----------------------------
  perturb_reps <- NULL
  SE_RMST <- NULL
  CI_RMST <- NULL
  p_value <- NULL
  z_value <- NULL

  if (isTRUE(inference)) {
    if (!is.null(seed)) set.seed(seed)

    n <- length(delta_i)
    perturb_reps <- numeric(B)

    for (b in seq_len(B)) {
      # Exponential perturbation weights with mean ~1 after normalization
      w <- stats::rexp(n, rate = 1)
      w <- w / mean(w)

      perturb_reps[b] <- sum(w * delta_i) / sum(w)
    }

    SE_RMST <- stats::sd(perturb_reps)

    alpha <- 1 - ci_level
    zcrit <- stats::qnorm(1 - alpha / 2)

    CI_RMST <- c(
      ATE - zcrit * SE_RMST,
      ATE + zcrit * SE_RMST
    )
    names(CI_RMST) <- c("lower", "upper")

    if (!is.na(SE_RMST) && SE_RMST > 0) {
      z_value <- ATE / SE_RMST
      p_value <- 2 * stats::pnorm(-abs(z_value))
    } else {
      z_value <- NA_real_
      p_value <- NA_real_
    }
  }

  res <- list(
    ATE_RMST = ATE,
    mean_RMST_Treated = mean(rmst_1),
    mean_RMST_Control = mean(rmst_0),
    tau = tau,
    patient_rmst_treated = rmst_1,
    patient_rmst_control = rmst_0,
    patient_delta_rmst = delta_i,
    inference = inference,
    B = if (isTRUE(inference)) B else NULL,
    SE_RMST = SE_RMST,
    CI_RMST = CI_RMST,
    z_value = z_value,
    p_value = p_value,
    # p_value = base::format.pval(p_value, digits = 3, eps = 1e-16),
    perturb_reps = perturb_reps
  )

  msg <- sprintf("Adjusted Delta RMST at tau = %s: %s time units",
                 tau, round(ATE, 3))

  if (isTRUE(inference) && !is.null(SE_RMST)) {
    msg <- paste0(
      msg,
      " | SE = ", round(SE_RMST, 3),
      " | ", round(ci_level * 100), "% CI = [",
      round(CI_RMST[1], 3), ", ", round(CI_RMST[2], 3), "]"
    )
  }

  message(msg)

  return(res)
}








#' Plot Adjusted Marginal RMST Contrast Over Time
#'
#' Generates a curve showing how the adjusted marginal restricted mean survival
#' time (RMST) contrast evolves across a sequence of restriction times.
#'
#' If \code{inference = TRUE}, the function additionally displays perturbation-based
#' Wald confidence intervals at each value of \code{tau}.
#'
#' @param fit A fitted \code{SuperSurv} ensemble object.
#' @param data A \code{data.frame} containing the covariates and the binary grouping variable.
#' @param trt_col Character string. The exact name of the binary grouping variable in \code{data}.
#' @param times Numeric vector of time points matching the prediction grid.
#' @param tau_seq Numeric vector. A sequence of restriction times (\code{tau}) to evaluate and plot.
#' @param inference Logical; if \code{TRUE}, compute perturbation-based confidence intervals.
#'   Defaults to \code{FALSE}.
#' @param B Integer. Number of perturbation replicates used when \code{inference = TRUE}.
#'   Defaults to \code{200}.
#' @param seed Optional integer seed for reproducibility.
#' @param ci_level Numeric scalar in \code{(0,1)} specifying the confidence level for the
#'   confidence interval. Defaults to \code{0.95}.
#'
#' @examples
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
#'   tau_seq = tau_grid,
#'   inference = TRUE,
#'   B = 100,
#'   seed = 123
#' )
#'
#' @return A \code{ggplot} object visualizing the adjusted marginal RMST contrast curve.
#' @export
plot_marginal_rmst_curve <- function(fit, data, trt_col, times, tau_seq,
                                     inference = FALSE, B = 200,
                                     seed = NULL, ci_level = 0.95) {
  requireNamespace("ggplot2", quietly = TRUE)

  results <- lapply(seq_along(tau_seq), function(i) {
    t <- tau_seq[i]

    # optional varying seed across tau values for reproducibility
    seed_i <- if (!is.null(seed)) seed + i - 1 else NULL

    res <- estimate_marginal_rmst(
      fit = fit,
      data = data,
      trt_col = trt_col,
      times = times,
      tau = t,
      inference = inference,
      B = B,
      seed = seed_i,
      ci_level = ci_level
    )

    out <- data.frame(
      Tau = t,
      Delta_RMST = if (!is.null(res$Delta_RMST)) res$Delta_RMST else res$ATE_RMST,
      stringsAsFactors = FALSE
    )

    if (isTRUE(inference)) {
      out$SE_RMST <- res$SE_RMST
      out$lower <- res$CI_RMST["lower"]
      out$upper <- res$CI_RMST["upper"]
      out$p_value <- res$p_value
    }

    out
  })

  res_df <- do.call(rbind, results)

  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = Tau, y = Delta_RMST))

  if (isTRUE(inference)) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.4,
        fill = "#e63946"
      )
  }

  p <- p +
    ggplot2::geom_line(color = "#e63946", linewidth = 1) +
    ggplot2::geom_point(color = "#1d3557", size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Adjusted Marginal RMST Contrast Over Time",
      subtitle = if (isTRUE(inference)) {
        paste0("Difference in RMST with ", round(ci_level * 100), "% perturbation-based confidence intervals")
      } else {
        "Difference in Restricted Mean Survival Time between groups"
      },
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
