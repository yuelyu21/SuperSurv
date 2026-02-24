#' Survival Probability Heatmap
#'
#' @param object A fitted SuperSurv object
#' @param newdata Test covariates (e.g., `X_te[1:50, ]`)
#' @param times The time grid to visualize
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @export
plot_survival_heatmap <- function(object, newdata, times) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # 1. Generate Predictions
  preds <- predict(object, newdata = newdata, new.times = times)$event.SL.predict

  # 2. Prepare Data: Sort patients by their final survival probability
  # This makes the heatmap look "organized" (High risk at top, Low risk at bottom)
  order_idx <- order(preds[, ncol(preds)])
  preds_sorted <- preds[order_idx, ]

  df_heatmap <- as.data.frame(preds_sorted)
  colnames(df_heatmap) <- times
  df_heatmap$PatientID <- factor(1:nrow(df_heatmap))

  # 3. Pivot to Long Format
  df_long <- df_heatmap %>%
    tidyr::pivot_longer(cols = -PatientID, names_to = "Time", values_to = "Survival") %>%
    dplyr::mutate(Time = as.numeric(Time))

  ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = PatientID, fill = Survival)) +
    ggplot2::geom_tile() +
    # Using the 'Plasma' palette to match the purple/red/yellow theme of survex
    ggplot2::scale_fill_viridis_c(option = "plasma", direction = 1, name = "S(t)") +
    ggplot2::labs(
      title = "SuperSurv: Patient Risk Trajectories",
      subtitle = "Sorted by predicted survival at final follow-up",
      x = "Time (Months)",
      y = "Individual Patients"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      # Match the dark purple text from your importance plots
      text = ggplot2::element_text(color = "#371ea3"),
      plot.title = ggplot2::element_text(face = "bold", size = 18)
    )
}
