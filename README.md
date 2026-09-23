# SuperSurv: A Unified Framework for Machine Learning Ensembles in Survival Analysis <img src="man/figures/logo.png" align="right" width="140"/>

[![R-CMD-check](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**SuperSurv** is an R package for building, evaluating, and interpreting ensemble models for right-censored survival data.

At its core, the package implements a Super Learner-style ensemble framework for continuous-time survival prediction under right censoring. Using inverse probability of censoring weighting (IPCW), it combines heterogeneous base learners by minimizing cross-validated prediction risk. The framework supports learners that return full survival curves as well as learners that return only risk scores, which are calibrated to a common survival-probability scale on a shared evaluation time grid.

Beyond ensemble fitting, **SuperSurv** provides tools for:

- hyperparameter tuning
- high-dimensional screening
- time-dependent model benchmarking
- SHAP-based interpretability
- covariate-adjusted restricted mean survival time (RMST) contrasts

The package also provides a more user-friendly model interface through `print()`, `summary()`, `coef()`, and exported accessors such as `event_weights()`, `censor_weights()`, `learner_names()`, `training_variables()`, `selected_variables()`, and `eval_times()`.

---

## 🚀 Get Started

The best place to start is the installation and setup tutorial:

👉 **[Tutorial 0: Installation & Setup](https://yuelyu21.github.io/SuperSurv/articles/installation.html)**

You can also browse the full documentation site here:

👉 **[SuperSurv website](https://yuelyu21.github.io/SuperSurv/)**

---

## 📦 Installation

Install the CRAN release:

```r
install.packages("SuperSurv")
```

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("yuelyu21/SuperSurv")
```

---

## 📦 Dependency philosophy

To keep installation lightweight, several heavier machine-learning engines are listed in `Suggests` rather than imported as strict dependencies. This means users can install only the modeling backends they plan to use. If a requested learner is unavailable, **SuperSurv** will prompt the user to install the required package.

---

## 📚 Included learners and screeners

**SuperSurv** currently standardizes a broad set of prediction wrappers and screening methods within a unified interface.

### Prediction learners

- **Machine learning:** random forests, generalized random forests, gradient boosting, XGBoost, support vector machines, BART, ranger, DeepSurv, DeepHit
- **Penalized/high-dimensional:** elastic net, ridge regression, CoxBoost, component-wise Cox boosting
- **Tree-based:** RPART
- **Classical/parametric:** Cox proportional hazards, Weibull, exponential, log-logistic, log-normal, generic parametric models, and generalized gamma/Gompertz/gamma through `surv.flexsurvreg`
- **Smoothing/splines:** generalized additive models, flexible parametric spline models, and penalized smooth hazard models through `surv.survPen`
- **Baseline models:** Kaplan–Meier

### Screening methods

- keep all features
- marginal Cox screening
- variance-based screening
- elastic-net-based screening
- random forest variable hunting

The framework is extensible, and users can add custom learners and screeners. See the extensibility vignette for details.

DeepSurv, DeepHit, and Cox-Time are experimental optional adapters. They require
`survivalmodels`, `reticulate`, and a user-managed Python environment containing
`torch`, `torchtuples`, and `pycox`; they currently accept only uniform
observation weights. They support prediction from fitted objects in the active
R/Python session; plain `saveRDS()` does not preserve their Python model objects.
`surv.survPen` also accepts only uniform observation weights. Its one-sided
formula interface allows nonlinear and time-varying effects; `surv.flexsurvreg`
supports native observation weights. All of these backends remain optional.

---

## 📖 Documentation

The package website includes tutorials covering:

- **00. Installation & Setup**
- **01. The SuperSurv Ensemble**
- **02. Model Performance**
- **03. Selection vs. Ensemble**
- **04. Screening Methods**
- **05. Hyperparameter Tuning**
- **06. Random Forests**
- **07. Parametric Models**
- **08. SHAP Interpretability**
- **09. Causal Inference (RMST)**
- **10. Extending SuperSurv**

---

## 📖 Citation

To cite the package, use:

```r
citation("SuperSurv")
```

If you would also like to cite the accompanying preprint:

Lyu, Y., Lin, S. H., Huang, X., & Li, Z. (2026).  
*SuperSurv: A Unified Framework for Machine Learning Ensembles in Survival Analysis*.  
bioRxiv.  
https://doi.org/10.64898/2026.03.11.711010

Related methodological work:

Westling, T., Luedtke, A., Gilbert, P. B., & Carone, M. (2024).  
*Inference for treatment-specific survival curves using machine learning*.  
Journal of the American Statistical Association.  
https://doi.org/10.1080/01621459.2023.2205060
