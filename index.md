# SuperSurv: Unified Ensemble Machine Learning for Survival Data ![](reference/figures/logo.png)

[![R-CMD-check](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml)
\## Overview

**SuperSurv** provides a unified, explainable machine learning ecosystem
for right-censored time-to-event data.

At its core, `SuperSurv` implements the **Double Super Learner**
framework. By utilizing Inverse Probability of Censoring Weighting
(IPCW) generated from a censoring-specific meta-learner, it constructs
mathematically rigorous, optimal convex combinations of parametric,
semi-parametric, and machine learning base algorithms (e.g., Random
Forests, XGBoost, Elastic Net) to minimize cross-validated risk.

Beyond the core ensemble engine, `SuperSurv` acts as a complete
Explainable AI (XAI) pipeline for survival analysis, featuring: \*
**High-Dimensional Variable Screening:** Automated Marginal, Variance,
and Elastic Net screening for genomic and high-dimensional clinical
data. \* **Automated Hyperparameter Grids:** Dynamic tuning for complex
algorithms. \* **Rigorous Performance Evaluation:** Native computation
of Time-Dependent Brier Scores, Uno’s C-index, and clinical calibration
metrics. \* **State-of-the-Art Interpretability:** Seamless integration
with Kernel SHAP and time-dependent explanations to demystify complex
survival predictions.

## Documentation & Tutorials

**[Visit the official SuperSurv
website](https://yuelyu21.github.io/SuperSurv)** for comprehensive
documentation, function references, and a complete suite of 8 in-depth
tutorials, including: 1. The SuperSurv Ensemble 2. Selection
vs. Ensemble Methods 3. Model Performance & Calibration 4. Screening
Methods 5. SHAP Interpretability 6. Scaling Up with Parallel Processing

## Installation

You can install the development version of `SuperSurv` directly from
GitHub using:

\`\`\`r \# install.packages(“devtools”)
devtools::install_github(“yuelyu21/SuperSurv”)

library(SuperSurv)

# Load built-in data

data(“metabric”, package = “SuperSurv”)

# Define predictor matrix and time grid

X \<- metabric\[, grep(“^x”, names(metabric))\] new.times \<- seq(50,
200, by = 25)

# Fit the Double Super Learner

fit \<- SuperSurv( time = metabric$duration,event = metabric$event, X =
X, newX = X, new.times = new.times, event.SL.library = c(“surv.coxph”,
“surv.weibull”, “surv.rfsrc”), cens.SL.library = c(“surv.coxph”),
parallel = TRUE, nFolds = 5 )
