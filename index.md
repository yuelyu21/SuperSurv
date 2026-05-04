# SuperSurv: A Unified Framework for Machine Learning Ensembles in Survival Analysis

[![R-CMD-check](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yuelyu21/SuperSurv/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**SuperSurv** is an R package for building, evaluating, and interpreting
ensemble models for right-censored survival data.

At its core, the package implements a Super Learner-style ensemble
framework for continuous-time survival prediction under right censoring.
Using inverse probability of censoring weighting (IPCW), it combines
heterogeneous base learners by minimizing cross-validated prediction
risk. The framework supports learners that return full survival curves
as well as learners that return only risk scores, which are calibrated
to a common survival-probability scale on a shared evaluation time grid.

Beyond ensemble fitting, **SuperSurv** provides tools for:

- hyperparameter tuning
- high-dimensional screening
- time-dependent model benchmarking
- SHAP-based interpretability
- covariate-adjusted restricted mean survival time (RMST) contrasts

The package also provides a more user-friendly model interface through
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), and exported accessors
such as
[`event_weights()`](https://yuelyu21.github.io/SuperSurv/reference/event_weights.md),
[`censor_weights()`](https://yuelyu21.github.io/SuperSurv/reference/event_weights.md),
[`learner_names()`](https://yuelyu21.github.io/SuperSurv/reference/learner_names.md),
[`training_variables()`](https://yuelyu21.github.io/SuperSurv/reference/training_variables.md),
[`selected_variables()`](https://yuelyu21.github.io/SuperSurv/reference/selected_variables.md),
and
[`eval_times()`](https://yuelyu21.github.io/SuperSurv/reference/eval_times.md).

------------------------------------------------------------------------

## 🚀 Get Started

The best place to start is the installation and setup tutorial:

👉 **[Tutorial 0: Installation &
Setup](https://yuelyu21.github.io/SuperSurv/articles/installation.html)**

You can also browse the full documentation site here:

👉 **[SuperSurv website](https://yuelyu21.github.io/SuperSurv/)**

------------------------------------------------------------------------

## 📦 Installation

Install the CRAN release:

``` r

install.packages("SuperSurv")
```

Install the development version from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("yuelyu21/SuperSurv")
```

------------------------------------------------------------------------

## 📦 Dependency philosophy

To keep installation lightweight, several heavier machine-learning
engines are listed in `Suggests` rather than imported as strict
dependencies. This means users can install only the modeling backends
they plan to use. If a requested learner is unavailable, **SuperSurv**
will prompt the user to install the required package.

------------------------------------------------------------------------

## 📚 Included learners and screeners

**SuperSurv** currently standardizes a broad set of prediction wrappers
and screening methods within a unified interface.

### Prediction learners

- **Machine learning:** random forests, gradient boosting, XGBoost,
  support vector machines, BART, ranger
- **Penalized/high-dimensional:** elastic net, ridge regression,
  CoxBoost
- **Tree-based:** RPART
- **Classical/parametric:** Cox proportional hazards, Weibull,
  exponential, log-logistic, log-normal, generic parametric models
- **Smoothing/splines:** generalized additive models
- **Baseline models:** Kaplan–Meier

### Screening methods

- keep all features
- marginal Cox screening
- variance-based screening
- elastic-net-based screening
- random forest variable hunting

The framework is extensible, and users can add custom learners and
screeners. See the extensibility vignette for details.

------------------------------------------------------------------------

## 📖 Documentation

The package website includes tutorials covering:

- **0. Installation & Setup**
- **1. The SuperSurv Ensemble**
- **2. Model Performance**
- **3. Selection vs. Ensemble**
- **4. Screening Methods**
- **5. Hyperparameter Tuning**
- **6. Random Forests**
- **7. Parametric Models**
- **8. SHAP Interpretability**
- **9. Causal Inference (RMST)**
- **10. Parallel Processing**
- **11. Extending SuperSurv**

------------------------------------------------------------------------

## 📖 Citation

To cite the package, use:

``` r

citation("SuperSurv")
```

If you would also like to cite the accompanying preprint:

Lyu, Y., Lin, S. H., Huang, X., & Li, Z. (2026).  
*SuperSurv: A Unified Framework for Machine Learning Ensembles in
Survival Analysis*.  
bioRxiv.  
<https://doi.org/10.64898/2026.03.11.711010>

Related methodological work:

Westling, T., Luedtke, A., Gilbert, P. B., & Carone, M. (2024).  
*Inference for treatment-specific survival curves using machine
learning*.  
Journal of the American Statistical Association.  
<https://doi.org/10.1080/01621459.2023.2205060>
