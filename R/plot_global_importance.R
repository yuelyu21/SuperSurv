#' Plot Global Feature Importance for SuperSurv
#' @param title Plot title.
#' @param shap_values The output from explain_shap.SuperSurv
#' @param top_n Number of features to show (default 10)
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @export
plot_global_importance <- function(shap_values, title = "SuperSurv: Ensemble Feature Importance", top_n = 10) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # 1. Calculate Mean |SHAP|
  # We use drop=FALSE and as.matrix to ensure it works with different data types
  importance <- colMeans(abs(as.matrix(shap_values)))

  # 2. Organize for ggplot
  df_imp <- data.frame(
    Feature = names(importance),
    Importance = as.numeric(importance)
  ) %>%
    dplyr::arrange(desc(Importance)) %>%
    dplyr::slice(1:top_n) # Take top N features

  # 3. Plot
  ggplot2::ggplot(df_imp, ggplot2::aes(x = stats::reorder(Feature, Importance), y = Importance)) +
    ggplot2::geom_col(fill = "steelblue", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = paste("Based on Weighted Ensemble SHAP Values"),
      x = NULL,
      y = "Mean |SHAP| (Impact on Mortality Risk)"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
