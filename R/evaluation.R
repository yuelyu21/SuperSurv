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
#'   S_mat = fit$pred,
#'   times = times
#' )
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
#'   S_mat = fit$pred,
#'   times = times,
#'   eval_time = 100,
#'   method = "uno"
#' )
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
#'    S_mat = fit$pred,
#'    times = times
#'  )
#' @export
eval_timeROC <- function(time, event, S_mat, times) {

  time <- as.numeric(time)
  event <- as.numeric(event)
  S_mat <- as.matrix(S_mat)
  times <- as.numeric(times)

  n <- length(time)
  if (nrow(S_mat) != n) {
    stop("Number of rows in S_mat must equal length(time).")
  }
  if (ncol(S_mat) != length(times)) {
    stop("Number of columns in S_mat must equal length(times).")
  }

  # Clamp probabilities
  eps_prob <- 1e-6
  S_mat[S_mat < eps_prob] <- eps_prob
  S_mat[S_mat > 1 - eps_prob] <- 1 - eps_prob

  # Estimate censoring survival G(t) using KM for censoring
  cens <- 1 - event
  fit_G <- survival::survfit(survival::Surv(time, cens) ~ 1)
  t_fit <- fit_G$time
  s_fit <- fit_G$surv
  eps_G <- 1e-6

  Gfun <- function(t) {
    idx <- findInterval(t, t_fit, rightmost.closed = TRUE)
    out <- ifelse(idx == 0, 1, s_fit[idx])
    pmax(out, eps_G)
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

    w_case[is_case] <- 1 / Gfun(time[is_case])
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
#' @return An object of class \code{"SuperSurv_eval"} containing benchmark
#'   metrics for the ensemble and base learners.
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

  if (isTRUE(verbose)) {
    message("Generating predictions on test data...")
  }
  preds <- predict(object, newdata = newdata, new.times = eval_times)

  k_models <- dim(preds$event.library.predict)[3]
  model_names <- dimnames(preds$event.library.predict)[[3]]

  if (is.null(model_names)) {
    model_names <- paste0("Base_Learner_", seq_len(k_models))
  }

  all_names <- c("SuperSurv_Ensemble", model_names)

  results <- data.frame(
    Model = all_names,
    IBS = numeric(length(all_names)),
    Uno_C = numeric(length(all_names)),
    iAUC = numeric(length(all_names)),
    stringsAsFactors = FALSE
  )

  process_matrix <- function(S_mat) {
    ibs_res <- eval_brier(time, event, S_mat, eval_times)$ibs

    c_res <- eval_cindex(
      time,
      event,
      S_mat,
      eval_times,
      eval_time = risk_time,
      method = "uno"
    )

    iauc_res <- suppressMessages(
      eval_timeROC(time, event, S_mat, eval_times)$iAUC
    )

    c(ibs_res, c_res, iauc_res)
  }

  if (isTRUE(verbose)) {
    message("Evaluating SuperSurv ensemble...")
  }
  results[1, 2:4] <- process_matrix(preds$event.predict)

  for (i in seq_len(k_models)) {
    if (isTRUE(verbose)) {
      message(sprintf(
        "Evaluating base learner %d/%d: %s",
        i, k_models, model_names[i]
      ))
    }
    results[i + 1, 2:4] <- process_matrix(preds$event.library.predict[, , i])
  }

  results$IBS <- round(results$IBS, 4)
  results$Uno_C <- round(results$Uno_C, 4)
  results$iAUC <- round(results$iAUC, 4)

  attr(results, "eval_times") <- eval_times
  attr(results, "risk_time") <- risk_time
  class(results) <- c("SuperSurv_eval", class(results))

  results
}








