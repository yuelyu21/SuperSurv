# 1. Appease CRAN checks for unquoted variables in ggplot2 and dplyr
utils::globalVariables(c(
  "Brier", "CD_AUC", "Direction", "Feature", "FeatureValue",
  "Feature_Value", "Importance", "Model", "PatientID", "SHAP",
  "SHAP_Value", "Survival", "Time", "mean_abs", "var"
))

# 2. Explicitly import base R stats functions and dplyr/magrittr operators
#' @importFrom stats as.formula coef median predict quantile var
#' @importFrom dplyr desc
#' @importFrom magrittr %>%
NULL
