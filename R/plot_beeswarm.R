#' Beeswarm Summary Plot for SuperSurv SHAP
#' 
#' @param shap_values The output from explain_shap.SuperSurv
#' @param data The covariate data used (X_explain)
#' @param top_n Number of features to display
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @export
plot_beeswarm <- function(shap_values, data, top_n = 10) {
  
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("ggforce", quietly = TRUE)
  
  # 1. Convert to standard data frame to avoid tibble issues
  shap_df <- as.data.frame(shap_values)
  data_df <- as.data.frame(data)
  
  # 2. Calculate overall importance to rank the Y-axis
  # We do this BEFORE pivoting to find the top features
  importance_df <- data.frame(
    Feature = names(shap_df),
    mean_abs = colMeans(abs(shap_df))
  ) %>%
    dplyr::arrange(desc(mean_abs)) %>%
    dplyr::slice(1:min(top_n, ncol(shap_df))) # EXPLICIT CALL
  
  top_features <- importance_df$Feature
  
  # 3. Reshape SHAP values to long format (only for top features)
  shap_long <- shap_df[, top_features, drop = FALSE] %>%
    tidyr::pivot_longer(cols = dplyr::everything(), 
                        names_to = "Feature", 
                        values_to = "SHAP")
  
  # 4. Reshape original data and scale (0 to 1) for color coding
  # We use a small epsilon (1e-9) to avoid division by zero
  data_long <- data_df[, top_features, drop = FALSE] %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), 
                                ~ (. - min(., na.rm=T)) / (max(., na.rm=T) - min(., na.rm=T) + 1e-9))) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), 
                        names_to = "Feature", 
                        values_to = "FeatureValue")
  
  # 5. Merge
  df_plot <- data.frame(
    Feature = shap_long$Feature,
    SHAP = shap_long$SHAP,
    FeatureValue = data_long$FeatureValue
  )
  
  # Ensure the Y-axis follows the importance ranking
  df_plot$Feature <- factor(df_plot$Feature, levels = rev(top_features))
  
  # 6. Plot
  ggplot2::ggplot(df_plot, ggplot2::aes(x = SHAP, y = Feature, color = FeatureValue)) +
    ggforce::geom_sina(alpha = 0.7, size = 1.8, method = "counts", maxwidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linetype = "solid", alpha = 0.5) +
    # SHAP standard colors: Blue (Low) to Red (High)
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#FF0D57", 
                                  breaks = c(0, 1), labels = c("Low", "High"),
                                  name = "Feature Value") +
    ggplot2::labs(
      title = "SuperSurv SHAP Summary",
      subtitle = paste("Top", length(top_features), "Features by Mean Absolute SHAP"),
      x = "SHAP Value (Impact on Mortality Risk)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold")
    )
}