#' Plot Comprehensive Survival Benchmarks
#' 
#' Generates beautiful ggplot2 visualizations for time-dependent Brier scores 
#' and Time-Dependent AUC without relying on external explainer packages.
#'
#' @param object A fitted SuperSurv object.
#' @param newdata Test covariates.
#' @param time True event times for the test set.
#' @param event True event indicators for the test set.
#' @param eval_times Numeric vector of evaluation times.
#' @return A list containing three elements:
#' \itemize{
#'   \item \code{brier_plot}: A \code{ggplot} object showing the IPCW Brier score over time.
#'   \item \code{auc_plot}: A \code{ggplot} object showing the Time-Dependent AUC.
#'   \item \code{plot_data}: A \code{data.frame} containing the raw metrics used for plotting.
#' }
#' @export
plot_benchmark <- function(object, newdata, time, event, eval_times) {
  
  requireNamespace("ggplot2", quietly = TRUE)
  
  # Official DALEX/survex color palette you love
  survex_colors <- c("#371ea3", "#f05a71", "#8bdcbe", "#ae2c87", "#ffa58c", "#4378bf", "#cccccc")
  
  # 1. Get Predictions natively
  preds <- predict(object, newdata = newdata, new.times = eval_times)
  
  k_models <- dim(preds$event.library.predict)[3]
  model_names <- dimnames(preds$event.library.predict)[[3]]
  all_names <- c("SuperSurv_Ensemble", model_names)
  
  plot_df <- data.frame()
  
  # 2. Loop through and compute native math for the plots
  for (i in seq_along(all_names)) {
    name <- all_names[i]
    message(paste("Calculating plotting data for:", name))
    
    # Grab the right matrix
    if (name == "SuperSurv_Ensemble") {
      S_mat <- preds$event.SL.predict
    } else {
      S_mat <- preds$event.library.predict[,, i - 1]
    }
    
    # Run our native math functions!
    brier_res <- eval_brier(time, event, S_mat, eval_times)
    roc_res <- suppressMessages(eval_timeROC(time, event, S_mat, eval_times))
    
    # Bind to the plotting dataframe
    temp_df <- data.frame(
      Time   = eval_times,
      Model  = name,
      Brier  = brier_res$brier_scores,
      CD_AUC = roc_res$AUC_curve
    )
    plot_df <- rbind(plot_df, temp_df)
  }
  
  # 3. Time-dependent Line Plots (Your exact ggplot code!)
  p_brier_td <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = Brier, color = Model)) +
    ggplot2::geom_line(linewidth = 1.2) + ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = survex_colors) +
    ggplot2::theme_minimal() + 
    ggplot2::labs(title = "IPCW Brier Score Over Time", y = "Brier Score (Lower is better)")
  
  p_auc_td   <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = CD_AUC, color = Model)) +
    ggplot2::geom_line(linewidth = 1.2) + ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = survex_colors) +
    ggplot2::theme_minimal() + 
    ggplot2::labs(title = "Time-Dependent AUC", y = "AUC (Higher is better)")
  
  return(list(
    brier_plot = p_brier_td, 
    auc_plot = p_auc_td,
    plot_data = plot_df # Returning the data just in case the user wants to customize!
  ))
}





#' Plot Predicted Survival Curves
#'
#' Generates a step-function plot of the predicted survival probabilities 
#' over time for a specific patient.
#'
#' @param preds A list containing predictions, specifically the object returned 
#'   by the \code{predict.SuperSurv} function.
#' @param patient_idx Integer. The row index of the patient in the test dataset 
#'   whose survival curve you want to plot. Defaults to 1.
#' @param ... Additional graphical parameters passed to the base \code{plot} 
#'   function (e.g., \code{col}, \code{lwd}, \code{lty}, \code{main}).
#' 
#' @return Called for its side effect of drawing a plot on the current graphics device. 
#'   Returns \code{NULL} invisibly.
#' @export
plot_predict <- function(preds, patient_idx = 1, ...) {
  
  # Extract times and the specific patient's survival probabilities
  times <- as.numeric(colnames(preds$event.SL.predict))
  surv_probs <- preds$event.SL.predict[patient_idx, ]
  
  # Base R step plot
  plot(times, surv_probs, type = "s", ylim = c(0, 1), 
       xlab = "Time", ylab = "Predicted Survival Probability",
       main = paste("Predicted Survival for Patient", patient_idx),
       col = "blue", lwd = 2, ...)
  
  invisible(NULL)
}





#' Plot IPCW Brier Score Curve
#'
#' Generates a line plot showing how the Inverse Probability of Censoring 
#' Weighted (IPCW) Brier score changes over the evaluation time grid.
#'
#' @param brier_eval_object A list containing Brier evaluation metrics, 
#'   specifically the object returned by the \code{eval_brier} function.
#' @param ... Additional graphical parameters passed to the base \code{plot} 
#'   function (e.g., \code{col}, \code{lwd}, \code{lty}, \code{main}).
#' 
#' @return Called for its side effect of drawing a plot on the current graphics device. 
#'   Returns \code{NULL} invisibly.
#' @export
plot_brier_curve <- function(brier_eval_object, ...) {
  
  # Extract times and scores from the evaluation object
  times <- brier_eval_object$times
  scores <- brier_eval_object$brier_scores
  
  # Base R line plot
  plot(times, scores, type = "l", ylim = c(0, max(scores, 0.3, na.rm = TRUE)),
       xlab = "Time", ylab = "IPCW Brier Score",
       main = "Brier Score over Time (Lower is Better)",
       col = "red", lwd = 2, ...)
  
  invisible(NULL)
}
