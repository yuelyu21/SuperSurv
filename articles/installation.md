# 0. Installation & Setup

## Welcome to SuperSurv!

`SuperSurv` is designed to be a unified ecosystem for machine learning
and survival analysis. However, installing 19 different machine learning
engines at once can take a long time and cause dependency conflicts on
some operating systems.

To make your experience as smooth as possible, `SuperSurv` uses a
**Modular Dependency Philosophy**.

The core package is incredibly lightweight and installs in seconds.
Heavy machine learning libraries (like XGBoost or Elastic Net) are only
required when you explicitly ask to use them!

------------------------------------------------------------------------

## Step 1: Install the Core Package

You can install the development version of `SuperSurv` directly from
GitHub using the `devtools` or `remotes` package:

``` r
# Install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install SuperSurv from GitHub
devtools::install_github("yuelyu21/SuperSurv")
```

Once installed, you can load the package and view all the available
modeling and screening wrappers:

``` r
library(SuperSurv)

# See all 19 prediction models and 6 screening algorithms!
list_wrappers()
```

------------------------------------------------------------------------

## Step 2: Install Base Learners (Optional but Recommended)

`SuperSurv` relies on external CRAN packages to run its various base
learners. If you try to run `surv.rfsrc` without having the
`randomForestSRC` package installed, `SuperSurv` will gently pause and
remind you to install it.

If you want to unlock the full power of the package right now, you can
copy and paste the following script to install the most commonly used
machine learning and interpretability engines:

``` r
# List of highly recommended modeling engines
ml_packages <- c(
  "survival",        # Classical Cox models
  "randomForestSRC", # Random Survival Forests
  "ranger",          # Fast Random Forests
  "xgboost",         # Extreme Gradient Boosting
  "glmnet",          # Elastic Net & Penalized Regression
  "rpart",           # Decision Trees
  "survex",          # Time-Dependent XAI (Interpretability)
  "fastshap"         # Kernel SHAP support
)

# Identify which ones you are missing
missing_pkgs <- ml_packages[!(ml_packages %in% installed.packages()[,"Package"])]

# Install the missing ones
if(length(missing_pkgs)) install.packages(missing_pkgs)
```

### Specialized Packages

A few wrappers require specialized packages that you might only need for
niche use cases: \* `surv.svm`: Requires `survivalsvm` \* `surv.gam`:
Requires `mgcv` \* `surv.coxboost`: Requires `CoxBoost`

------------------------------------------------------------------------

## Step 3: You’re Ready!

Your environment is now completely set up. You are ready to build your
first optimal survival ensemble!

👉 **[Click here to proceed to Tutorial 1: The SuperSurv
Ensemble](https://yuelyu21.github.io/SuperSurv/articles/supersurv-ensemble.md)**
