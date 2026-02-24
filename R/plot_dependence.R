#' Plot SHAP Dependence for SuperSurv
#' 
#' @param shap_values The output from explain_shap.SuperSurv
#' @param data The original covariate data used for the explanation (X_explain)
#' @param feature_name String name of the column to plot
#' @param title Optional custom title.
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @export
plot_dependence <- function(shap_values, data, feature_name, title = NULL) {
  
  requireNamespace("ggplot2", quietly = TRUE)
  if(is.null(title)) title <- paste("SHAP Dependence:", feature_name)
  
  df_plot <- data.frame(
    Feature_Value = data[[feature_name]],
    SHAP_Value = shap_values[[feature_name]]
  )
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = Feature_Value, y = SHAP_Value)) +
    # 1. Smooth line in background
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "black", 
                         se = TRUE, fill = "gray80", linetype = "dashed") +
    # 2. Points on top
    ggplot2::geom_point(alpha = 0.7, color = "firebrick", size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    ggplot2::labs(title = title, x = paste(feature_name, "Value"), y = "SHAP Value") +
    ggplot2::theme_minimal(base_size = 14)
}