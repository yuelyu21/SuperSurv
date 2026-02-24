#' Waterfall Plot for an Individual Patient
#' 
#' @param shap_values The output from explain_shap.*
#' @param patient_index The row index of the patient to explain
#' @param top_n Number of features to show (default 10)
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @export
plot_patient_waterfall <- function(shap_values, patient_index = 1, top_n = 10) {
  
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  # 1. Extract the SHAP values for that specific patient
  row_vals <- as.data.frame(shap_values)[patient_index, ]
  
  df_p <- data.frame(
    Feature = names(row_vals),
    SHAP = as.numeric(row_vals)
  ) %>%
    dplyr::arrange(desc(abs(SHAP))) %>%
    dplyr::slice(1:top_n) %>%
    dplyr::mutate(Direction = ifelse(SHAP > 0, "Higher Risk", "Lower Risk"))
  
  # 2. Plot
  ggplot2::ggplot(df_p, ggplot2::aes(x = stats::reorder(Feature, SHAP), y = SHAP, fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Higher Risk" = "#FF0D57", "Lower Risk" = "#1E88E5")) +
    ggplot2::labs(
      title = paste("Patient", patient_index, ": Risk Factor Contribution"),
      subtitle = "How individual features pushed this patient away from the average risk",
      x = NULL,
      y = "SHAP Value (Impact on Risk Score)"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}