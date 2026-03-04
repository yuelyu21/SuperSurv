# SuperSurv: A Unified Ecosystem for Machine Learning, Ensembles, and Interpretability in Survival Data

[![R-CMD-check](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**SuperSurv** provides a mathematically rigorous, unified machine learning ecosystem for right-censored time-to-event data. 

At its core, it implements an advanced **Ensemble Super Learner** framework. By utilizing Inverse Probability of Censoring Weighting (IPCW), it automatically constructs optimal convex combinations of parametric, semi-parametric, and machine learning base algorithms to minimize cross-validated risk.

Beyond the meta-learner, **SuperSurv** acts as a complete Explainable AI (XAI) pipeline, featuring automated hyperparameter tuning grids, rigorous time-dependent benchmarking, and seamless integration with XAI ecosystems like Kernel SHAP.

---

## 🚀 Let's Try the Package!

Ready to dive in? The best way to understand the power of SuperSurv is to see it in action. 

👉 **[Click here to get started with Tutorial 0: Installation & Setup!](https://yuelyu21.github.io/SuperSurv/articles/installation.html)**

---

## 📦 Dependencies Philosophy

To keep the core SuperSurv package lightweight and incredibly fast to install, heavy machine learning libraries (like XGBoost, Random Survival Forests, and Elastic Net) are listed as **Suggests** rather than strict requirements. 

This modular design means you only need to install the specific mathematical engines you actually plan to use! If you try to call a base learner that you haven't installed yet, SuperSurv will gently pause and remind you to install it. 

---

## 📚 The SuperSurv Library

SuperSurv standardizes the modeling API for **19 prediction algorithms** and **6 high-dimensional screening algorithms**. 

### Prediction Models
* **Machine Learning:** Random Forest, XGBoost, Support Vector Machines, Gradient Boosting, BART, Ranger
* **Penalized/High-Dimensional:** Elastic Net, Ridge Regression
* **Tree-Based:** RPART, CoxBoost
* **Parametric/Classical:** Cox Proportional Hazards, Weibull, Exponential, Log-Logistic, Log-Normal, Generic Parametric
* **Smoothing/Splines:** Generalized Additive Models (GAM)
* **Baselines:** Kaplan-Meier

### Screening Algorithms
Used to automatically filter massive datasets (like genomic data) before fitting complex models:
* Keep All Features
* Marginal Cox Screening
* Variance-based Screening
* Penalized Screening (Elastic Net)
* Random Forest Variable Hunting

---

## 📊 Comprehensive Documentation

Visit our [official website](https://yuelyu21.github.io/SuperSurv/) for a complete suite of in-depth tutorials:
* **0. Installation & Setup:** Get your R environment ready.
* **1. The SuperSurv Ensemble:** Build and train your first meta-learner.
* **2. Model Performance:** Evaluate Time-Dependent Brier Scores and Uno's C-index.
* **3. Selection vs. Ensemble:** Compare evaluation approaches.
* **4. Screening Methods:** Handle high-dimensional genomic data.
* **5. Hyperparameter Tuning:** Automate algorithmic grid searches.
* **6. Random Forests:** A deep dive into machine learning wrappers.
* **7. Parametric Models:** Classical statistical approaches.
* **8. SHAP Interpretability:** Demystify black-box models with global and local XAI.
* **9. Causal Inference (RMST):** Evaluate treatment effects over time.
* **10. Parallel Processing:** Scale up your computations for massive datasets.

---

## Installation

You can install the development version of `SuperSurv` directly from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("yuelyu21/SuperSurv")

