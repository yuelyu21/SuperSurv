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
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{brier_scores}: A numeric vector of Brier scores at each time point.
#'   \item \code{ibs}: The Integrated Brier Score over the range `\code{tmin}, \code{tmax}`.
#'   \item \code{times}: The time grid used.
#' }
#' @export
eval_brier <- function(time, event, S_mat, times, tmin = min(times), tmax = max(times)) {

  n <- length(time)
  T_len <- length(times)
  S_mat <- as.matrix(S_mat)

  if (nrow(S_mat) != n || ncol(S_mat) != T_len) {
    stop("S_mat dimensions (", nrow(S_mat), "x", ncol(S_mat),
         ") do not match n (", n, ") and length of times (", T_len, ").")
  }

  # Clamp probabilities to avoid exact 0/1 floating issues
  epsProb <- 1e-6
  epsG <- 1e-6
  S_mat[S_mat < epsProb] <- epsProb
  S_mat[S_mat > 1 - epsProb] <- 1 - epsProb

  # Helper: KM Censoring Survival G(t)
  cens <- 1 - event
  fit_G <- survival::survfit(survival::Surv(time, cens) ~ 1)
  t_fit <- fit_G$time
  s_fit <- fit_G$surv

  Gfun <- function(t) {
    idx <- findInterval(t, t_fit, rightmost.closed = TRUE)
    out <- ifelse(idx == 0, 1, s_fit[idx])
    pmax(out, epsG)
  }

  # Calculate Brier Score at each time t
  bs <- numeric(T_len)
  G_Ti <- Gfun(time) # G evaluated at each observed time

  for (j in seq_len(T_len)) {
    t <- times[j]
    S_t <- S_mat[, j]

    Y <- as.integer(time > t)
    G_t <- Gfun(t)

    # Graf IPCW Weights
    w <- ifelse(time <= t & event == 1, 1 / G_Ti,
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
#' @export
eval_cindex <- function(time, event, S_mat, times, eval_time, method = "uno") {

  requireNamespace("survival", quietly = TRUE)

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
#' Evaluates the cumulative/dynamic Time-Dependent AUC and Integrated AUC (iAUC)
#' using the \code{timeROC} package with IPCW adjustment.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param S_mat A numeric matrix of predicted survival probabilities.
#' @param times Numeric vector of evaluation times matching the columns of \code{S_mat}.
#'
#' @return A list containing the \code{AUC} at each time point and the \code{iAUC}.
#' @export
eval_timeROC <- function(time, event, S_mat, times) {

  requireNamespace("timeROC", quietly = TRUE)

  # For cumulative/dynamic AUC, the risk score is evaluated at baseline
  # or tracking over time. We can use the mean predicted risk across the time grid,
  # or the risk at the maximum time. A standard approach for matrix outputs
  # is using the risk score at the final prediction time.

  risk_score <- 1 - S_mat[, ncol(S_mat)]

  # Calculate Time-Dependent ROC using IPCW (marginal Kaplan-Meier for censoring)
  roc_fit <- timeROC::timeROC(
    T = time,
    delta = event,
    marker = risk_score,
    cause = 1,
    weighting = "marginal",
    times = times,
    iid = FALSE # Set to TRUE if you want confidence intervals later
  )

  # Integrated AUC (using the timeROC built-in summary)
  # timeROC doesn't explicitly return iAUC in the object, but we can compute it
  # via the trapezoidal rule over the returned AUCs just like IBS!

  auc_vals <- roc_fit$AUC
  valid_idx <- which(!is.na(auc_vals))

  if (length(valid_idx) < 2) {
    iauc <- NA
  } else {
    t_sub <- times[valid_idx]
    a_sub <- auc_vals[valid_idx]
    diffs <- diff(t_sub)
    heights <- (a_sub[-1] + a_sub[-length(a_sub)]) / 2
    iauc <- sum(diffs * heights) / (max(t_sub) - min(t_sub))
  }

  return(list(
    AUC_curve = roc_fit$AUC,
    times = roc_fit$times,
    iAUC = iauc
  ))
}




#' Evaluate SuperSurv Predictions on Test Data
#'
#' Computes the Integrated Brier Score (IBS), Uno's C-index, and Integrated AUC (iAUC)
#' for the SuperSurv ensemble and all individual base learners.
#'
#' @param object A fitted \code{SuperSurv} object.
#' @param newdata A data.frame of test covariates.
#' @param time Numeric vector of observed follow-up times for the test set.
#' @param event Numeric vector of event indicators for the test set.
#' @param eval_times Numeric vector of times at which to evaluate survival predictions.
#' @param risk_time Numeric. The specific time horizon to use when extracting risk
#'   scores for Uno's C-index. Defaults to the median of \code{eval_times}.
#'
#' @return A data.frame containing the benchmark metrics for all models.
#' @export
eval_summary <- function(object, newdata, time, event, eval_times, risk_time = median(eval_times)) {

  # 1. Generate Predictions for the Ensemble and Base Learners
  cat("\nGenerating predictions on test data...\n")
  preds <- predict(object, newdata = newdata, new.times = eval_times)

  # 2. Setup Storage
  k_models <- dim(preds$event.library.predict)[3]
  model_names <- dimnames(preds$event.library.predict)[[3]]
  all_names <- c("SuperSurv_Ensemble", model_names)

  results <- data.frame(
    Model = all_names,
    IBS = numeric(length(all_names)),
    Uno_C = numeric(length(all_names)),
    iAUC = numeric(length(all_names)),
    stringsAsFactors = FALSE
  )

  # 3. Internal Helper to compute all metrics for a single matrix
  process_matrix <- function(S_mat) {
    # IBS
    ibs_res <- eval_brier(time, event, S_mat, eval_times)$ibs

    # Uno's C-index
    c_res <- eval_cindex(time, event, S_mat, eval_times, eval_time = risk_time, method = "uno")

    # iAUC (Suppressing timeROC text output for a clean console)
    iauc_res <- suppressMessages(eval_timeROC(time, event, S_mat, eval_times)$iAUC)

    return(c(ibs_res, c_res, iauc_res))
  }

  # 4. Evaluate Ensemble
  cat("Evaluating SuperSurv Ensemble...\n")
  results[1, 2:4] <- process_matrix(preds$event.SL.predict)

  # 5. Evaluate Base Learners
  for (i in seq_len(k_models)) {
    cat(sprintf("Evaluating Base Learner %d/%d: %s...\n", i, k_models, model_names[i]))
    results[i+1, 2:4] <- process_matrix(preds$event.library.predict[,,i])
  }

  # 6. Format and Print the Benchmark Table
  results$IBS <- round(results$IBS, 4)
  results$Uno_C <- round(results$Uno_C, 4)
  results$iAUC <- round(results$iAUC, 4)

  cat("\n========================================================\n")
  cat("             SuperSurv Evaluation Benchmark             \n")
  cat("========================================================\n")
  print(results, row.names = FALSE)
  cat("========================================================\n")
  cat(sprintf("* IBS & iAUC integrated over: [%.2f, %.2f]\n", min(eval_times), max(eval_times)))
  cat(sprintf("* Uno's C-index evaluated using risk at time: %.2f\n", risk_time))
  cat("Note: Lower IBS is better. Higher Uno_C and iAUC are better.\n")

  invisible(results)
}
