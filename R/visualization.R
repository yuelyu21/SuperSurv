# Define global variables for ggplot2 to satisfy R CMD check
utils::globalVariables(c(
  "C_Index", "Mean_Predicted", "Observed", "Survival_Prob",
  "Patient", "Feature", "Importance", "SHAP", "FeatureValue",
  "Direction", "Time", "Survival",
  # add these:
  "Tau", "ATE", "Observed_Time", "Predicted_RMST", "Event"
))

#' Plot Longitudinal Benchmark Metrics
#'
#' Generates time-dependent performance curves comparing the SuperSurv ensemble
#' against its base learners, or evaluates a single standalone learner.
#'
#' @param object A fitted SuperSurv object OR a fitted standalone learner.
#' @param newdata A data.frame of test covariates.
#' @param time Numeric vector of observed follow-up times for the test set.
#' @param event Numeric vector of event indicators for the test set.
#' @param eval_times Numeric vector of times at which to evaluate predictions.
#' @param metrics Character vector specifying which plots to return.
#'   Options: "brier", "auc", "cindex". Defaults to all three.
#'
#' @return A combined patchwork ggplot object, or a single ggplot if only one metric is selected.
#' @export
plot_benchmark <- function(object, newdata, time, event, eval_times,
                           metrics = c("brier", "auc", "cindex")) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("patchwork", quietly = TRUE)

  message("Generating predictions for benchmark plots...")
  preds <- predict(object, newdata = newdata, new.times = eval_times)

  get_metrics_df <- function(model_name, S_mat) {
    brier_vals <- eval_brier(time, event, S_mat, eval_times)$brier_scores
    auc_vals <- suppressMessages(eval_timeROC(time, event, S_mat, eval_times)$AUC_curve)

    cindex_vals <- vapply(eval_times, function(t) {
      eval_cindex(time, event, S_mat, eval_times, eval_time = t, method = "uno")
    }, numeric(1))

    data.frame(
      Time = eval_times,
      Model = model_name,
      Brier = brier_vals,
      CD_AUC = auc_vals,
      C_Index = cindex_vals,
      stringsAsFactors = FALSE
    )
  }

  message("Calculating time-dependent metrics...")

  res_list <- list()

  # SMART EXTRACTION LOGIC
  if (is.list(preds) && !is.null(preds$event.predict)) {
    # SuperSurv ensemble
    k_models <- dim(preds$event.library.predict)[3]

    model_names <- object$event.libraryNames
    if (is.null(model_names) || length(model_names) != k_models) {
      # fall back to dimnames if present
      dn <- dimnames(preds$event.library.predict)[[3]]
      if (!is.null(dn) && length(dn) == k_models) {
        model_names <- dn
      } else {
        model_names <- paste0("Base_Learner_", seq_len(k_models))
      }
    }

    all_names <- c("SuperSurv_Ensemble", model_names)

    res_list[[1]] <- get_metrics_df("SuperSurv_Ensemble", preds$event.predict)
    for (i in seq_len(k_models)) {
      res_list[[i + 1]] <- get_metrics_df(model_names[i], preds$event.library.predict[, , i])
    }

  } else if (is.matrix(preds)) {
    # Standalone learner
    model_name <- class(object)[1]
    if (is.null(model_name) || model_name %in% c("matrix")) model_name <- "Standalone_Model"

    all_names <- model_name
    res_list[[1]] <- get_metrics_df(model_name, preds)

  } else {
    stop("Unrecognized prediction format. Must be a SuperSurv object or a valid base learner with matrix predictions.")
  }

  plot_df <- do.call(rbind, res_list)
  plot_df$Model <- factor(plot_df$Model, levels = all_names)

  survex_colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                     "#ff7f00", "#ffff33", "#a65628", "#f781bf")
  survex_colors <- rep(survex_colors, length.out = length(all_names))
  names(survex_colors) <- all_names

  plots <- list()

  if ("brier" %in% metrics) {
    plots$brier <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = Brier, color = Model)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = survex_colors) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "IPCW Brier Score Over Time", y = "Brier Score (Lower is better)")
  }

  if ("auc" %in% metrics) {
    plots$auc <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = CD_AUC, color = Model)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = survex_colors) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Time-Dependent AUC", y = "AUC (Higher is better)")
  }

  if ("cindex" %in% metrics) {
    plots$cindex <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = C_Index, color = Model)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = survex_colors) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Uno's C-index Over Time", y = "C-index (Higher is better)")
  }

  if (length(plots) == 0) stop("No valid metrics selected.")
  if (length(plots) == 1) return(plots[[1]])

  patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom",
                   legend.direction = "vertical") &
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    color = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}


#' Plot Predicted Survival Curves
#'
#' Generates a ggplot2 step-function plot of the predicted survival probabilities
#' over time for one or more specific patients.
#'
#' @param preds A list containing predictions, specifically the object returned
#'   by the \code{predict.SuperSurv} function.
#' @param eval_times Numeric vector of times at which predictions were evaluated.
#' @param patient_idx Integer vector. The row indices of the patients in the test
#'   dataset whose survival curves you want to plot. Defaults to 1.
#'
#' @return A \code{ggplot} object showing the survival curves.
#' @export
plot_predict <- function(preds, eval_times, patient_idx = 1) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("patchwork", quietly = TRUE) # Add this!

  # Extract the prediction matrix
  S_mat <- preds$event.predict

  # Initialize an empty data frame to store all patients' data
  plot_df <- data.frame(Time = numeric(), Survival_Prob = numeric(), Patient = character(),
                        stringsAsFactors = FALSE)

  # Loop through the requested patients and format the data
  for (i in patient_idx) {
    temp_df <- data.frame(
      Time = eval_times,
      Survival_Prob = S_mat[i, ],
      Patient = paste("Patient", i)
    )
    plot_df <- rbind(plot_df, temp_df)
  }

  # Convert Patient to a factor so the legend stays in numeric order
  plot_df$Patient <- factor(plot_df$Patient, levels = paste("Patient", patient_idx))

  # Generate the ggplot using geom_step (the ggplot equivalent of type = "s")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = Survival_Prob, color = Patient)) +
    ggplot2::geom_step(linewidth = 1.2) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Predicted Patient Survival Curves",
      x = "Follow-up Time",
      y = "Predicted Survival Probability"
    ) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold"))

  return(p)
}




#' Plot Survival Calibration Curve
#'
#' @param object A fitted SuperSurv object OR a standalone base learner.
#' @param newdata A data.frame of test covariates.
#' @param time Numeric vector of observed follow-up times for the test set.
#' @param event Numeric vector of event indicators for the test set.
#' @param eval_time Numeric. A single time point at which to assess calibration.
#' @param bins Integer. Defaults to 5.
#' @return A ggplot object.
#' @export
plot_calibration <- function(object, newdata, time, event, eval_time, bins = 5) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  preds <- predict(object, newdata = newdata, new.times = eval_time)

  if (is.list(preds) && !is.null(preds$event.predict)) {
    S_pred <- preds$event.predict[, 1]
  } else if (is.matrix(preds)) {
    S_pred <- preds[, 1]
  } else {
    stop("Unrecognized model object or prediction output.")
  }

  df <- data.frame(time = time, event = event, pred = S_pred)

  # Tie-safe quantile binning
  qs <- stats::quantile(df$pred, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE)
  qs <- unique(qs)
  if (length(qs) <= 2) stop("Not enough unique prediction values to form calibration bins.")

  df$bin <- cut(df$pred, breaks = qs, include.lowest = TRUE, labels = FALSE)

  calib_data <- data.frame(Bin = integer(), Mean_Predicted = numeric(), Observed = numeric())

  for (i in sort(unique(df$bin))) {
    sub_df <- df[df$bin == i, , drop = FALSE]
    mean_pred <- mean(sub_df$pred, na.rm = TRUE)

    km <- survival::survfit(survival::Surv(time, event) ~ 1, data = sub_df)
    idx <- findInterval(eval_time, km$time)
    obs_surv <- ifelse(idx == 0, 1, km$surv[idx])

    calib_data <- rbind(calib_data, data.frame(
      Bin = i,
      Mean_Predicted = mean_pred,
      Observed = obs_surv
    ))
  }

  ggplot2::ggplot(calib_data, ggplot2::aes(x = Mean_Predicted, y = Observed)) +
    ggplot2::geom_point(size = 4, color = "#e41a1c") +
    ggplot2::geom_line(color = "#e41a1c", linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "black",
                         linetype = "dashed", linewidth = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = sprintf("Calibration Curve at Time = %s", eval_time),
      subtitle = "Red line = Model Calibration | Dashed line = Perfect Calibration",
      x = "Predicted Survival Probability",
      y = "Observed Survival Probability (Kaplan-Meier)"
    )
}






#' Plot Predicted Survival Curves
#'
#' @param preds A list containing SuperSurv predictions OR a raw prediction matrix.
#' @param eval_times Numeric vector of times at which predictions were evaluated.
#' @param patient_idx Integer vector. Defaults to 1.
#' @return A ggplot object.
#' @export
plot_predict <- function(preds, eval_times, patient_idx = 1) {
  requireNamespace("ggplot2", quietly = TRUE)

  if (is.list(preds) && !is.null(preds$event.predict)) {
    S_mat <- preds$event.predict
  } else if (is.matrix(preds)) {
    S_mat <- preds
  } else {
    stop("preds must be a prediction list from SuperSurv or a prediction matrix.")
  }

  plot_df <- data.frame(Time = numeric(), Survival_Prob = numeric(), Patient = character(),
                        stringsAsFactors = FALSE)

  for (i in patient_idx) {
    temp_df <- data.frame(
      Time = eval_times,
      Survival_Prob = S_mat[i, ],
      Patient = paste("Patient", i)
    )
    plot_df <- rbind(plot_df, temp_df)
  }

  plot_df$Patient <- factor(plot_df$Patient, levels = paste("Patient", patient_idx))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = Survival_Prob, color = Patient)) +
    ggplot2::geom_step(linewidth = 1.2) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Predicted Patient Survival Curves",
      x = "Follow-up Time",
      y = "Predicted Survival Probability"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}


