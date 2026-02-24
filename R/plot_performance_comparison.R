#' Benchmark Ensemble vs. Baseline
#' 
#' @param explainers A list of survex explainers (e.g., list(Ensemble=e1, Cox=e2))
plot_performance_comparison <- function(explainers) {
  
  # 1. Calculate Model Performance (Brier Score & C-Index over time)
  perf <- lapply(explainers, survex::model_performance)
  
  # 2. Plot Brier Score (Lower is Better)
  # This shows which model is more 'accurate' at every month/day
  p1 <- plot(perf[[1]], perf[[2]], metrics = "brier_score") + 
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Prediction Error (Brier Score) Over Time")
  
  # 3. Plot C-Index (Higher is Better)
  p2 <- plot(perf[[1]], perf[[2]], metrics = "c_index") + 
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Discrimination (C-Index) Over Time")
  
  return(list(p1, p2))
}