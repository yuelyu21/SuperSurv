## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-core, eval=FALSE-------------------------------------------------
# # Install devtools if you don't have it
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# 
# # Install SuperSurv from GitHub
# devtools::install_github("yuelyu21/SuperSurv")

## ----load-package, eval=FALSE-------------------------------------------------
# library(SuperSurv)
# 
# # See all 19 prediction models and 6 screening algorithms!
# list_wrappers()

## ----install-suggests, eval=FALSE---------------------------------------------
# # List of highly recommended modeling engines
# ml_packages <- c(
#   "survival",        # Classical Cox models
#   "randomForestSRC", # Random Survival Forests
#   "ranger",          # Fast Random Forests
#   "xgboost",         # Extreme Gradient Boosting
#   "glmnet",          # Elastic Net & Penalized Regression
#   "rpart",           # Decision Trees
#   "survex",          # Time-Dependent XAI (Interpretability)
#   "fastshap"         # Kernel SHAP support
# )
# 
# # Identify which ones you are missing
# missing_pkgs <- ml_packages[!(ml_packages %in% installed.packages()[,"Package"])]
# 
# # Install the missing ones
# if(length(missing_pkgs)) install.packages(missing_pkgs)

