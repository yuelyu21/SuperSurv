# 4. High-Dimensional Data & Variable Screening

## Introduction

In modern clinical research and bioinformatics, it is common to
encounter high-dimensional datasets where the number of predictors
(genes, biomarkers, clinical features) far exceeds the number of
observed events ($p \gg n$).

Feeding thousands of raw variables directly into complex machine
learning algorithms introduces two massive problems: 1. **Computational
Bottlenecks:** Algorithms like Random Survival Forests become
prohibitively slow. 2. **Overfitting:** Models will mathematically
memorize background noise instead of true biological signals, destroying
out-of-sample prediction accuracy.

`SuperSurv` elegantly handles this through **Variable Screening**.
Screening acts as a computational filter *before* the models are fitted.
It evaluates all predictors, drops the statistical noise, and only
passes the most highly associated features into your base learners.

## 1. Simulating High-Dimensional Data

To demonstrate how screening works, we will use our standard `metabric`
dataset, but we will artificially inject 100 columns of pure random
noise. This simulates a high-dimensional genomic dataset where most
measured features have absolutely no relationship to patient survival.

``` r
library(SuperSurv)
library(survival)

data("metabric", package = "SuperSurv")
set.seed(42)

# Create 100 columns of random noise (simulating irrelevant genes/biomarkers)
n_patients <- nrow(metabric)
noise_data <- matrix(rnorm(n_patients * 20), nrow = n_patients)
colnames(noise_data) <- paste0("noise_", 1:20)
# Combine the real data with the noise
high_dim_data <- cbind(metabric, noise_data)

# Standard Train/Test Split
train_idx <- sample(1:nrow(high_dim_data), 0.7 * nrow(high_dim_data))
train <- high_dim_data[train_idx, ]
test  <- high_dim_data[-train_idx, ]

X_tr <- train[, grep("^x|^noise", names(high_dim_data))]
X_te <- test[, grep("^x|^noise", names(high_dim_data))]
new.times <- seq(50, 200, by = 25)
```

## 2. Defining the Screening Algorithms

Just as we have a library for prediction algorithms, `SuperSurv`
utilizes a library for screening algorithms. You can pair any prediction
model with any screening function to create a streamlined pipeline.

`SuperSurv` includes a comprehensive suite of built-in screeners: \*
**`screen.marg`**: Marginal screening (keeps variables with strong
univariate associations). \* **`screen.glmnet`**: Penalized Lasso
screening (shrinks irrelevant feature coefficients to exact zero). \*
**`screen.elasticnet`**: Elastic Net screening (combines Lasso and Ridge
penalties; highly recommended for highly correlated biological features
like gene expression data). \* **`screen.rfsrc`**: Tree-based screening
(keeps features with high Random Forest variable importance). \*
**`screen.var`**: Unsupervised screening (drops features with near-zero
variance). \* **`screen.all`**: A baseline passthrough (keeps all
variables without screening).

Let’s build a library that tests different combinations of prediction
and screening:

``` r
# We pass a list of vectors. The first element is the prediction model, 
# and the second element is the screening algorithm.

my_screen_library <- list(
  c("surv.coxph", "screen.all"),         # Baseline: Cox model with ALL variables
  c("surv.coxph", "screen.marg"),        # Screen by marginal association, then fit Cox
  c("surv.weibull", "screen.elasticnet"),# Screen via Elastic Net, then fit Weibull
  c("surv.rpart", "screen.var")          # Drop zero-variance noise, then fit a Tree
)

# For the censoring library, we use an unscreened approach 
cens_library <- c("surv.coxph")
```

## 3. Training with Embedded Screening

When we pass this paired library into `SuperSurv`, the meta-learner
handles the cross-validation with absolute statistical rigor.

**Crucial Methodological Detail:** The screening step is performed
*inside* the cross-validation folds. If you screen the entire dataset
before splitting it into folds, you cause “data leakage,” leading to
overly optimistic benchmark metrics. `SuperSurv` protects you from this
automatically.

``` r
fit_highdim <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newX = X_te,
  new.times = new.times,
  event.SL.library = my_screen_library, # Pass our screened library
  cens.SL.library = cens_library,
  control = list(saveFitLibrary = TRUE),
  verbose = T,
  nFolds = 3
)
```

## 4. Evaluating the Screened Ensemble

Despite injecting 100 columns of pure noise, our ensemble still performs
exceptionally well because the screening algorithms successfully
filtered out the irrelevant variables before the base learners could
overfit to them.

We can evaluate this exactly as we did in previous tutorials using
[`eval_summary()`](https://yuelyu21.github.io/SuperSurv/reference/eval_summary.md)
and
[`plot_benchmark()`](https://yuelyu21.github.io/SuperSurv/reference/plot_benchmark.md).

``` r
# Check the integrated performance
screened_performance <- eval_summary(
  object = fit_highdim,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_times = new.times
)
#> 
#> Generating predictions on test data...
#> Evaluating SuperSurv Ensemble...
#> Evaluating Base Learner 1/4: surv.coxph_screen.all...
#> Evaluating Base Learner 2/4: surv.coxph_screen.marg...
#> Evaluating Base Learner 3/4: surv.weibull_screen.elasticnet...
#> Evaluating Base Learner 4/4: surv.rpart_screen.var...
#> 
#> ========================================================
#>              SuperSurv Evaluation Benchmark             
#> ========================================================
#>                           Model    IBS  Uno_C   iAUC
#>              SuperSurv_Ensemble 0.2014 0.6373 0.6779
#>           surv.coxph_screen.all 0.2009 0.6350 0.6794
#>          surv.coxph_screen.marg 0.2004 0.6270 0.6775
#>  surv.weibull_screen.elasticnet 0.2005 0.6324 0.6798
#>           surv.rpart_screen.var 0.2106 0.6104 0.6349
#> ========================================================
#> * IBS & iAUC integrated over: [50.00, 200.00]
#> * Uno's C-index evaluated using risk at time: 125.00
#> Note: Lower IBS is better. Higher Uno_C and iAUC are better.

# Plot the benchmark
# plot_benchmark(fit_highdim, newdata = X_te, time = test$duration, event = test$event, eval_times = new.times)
```

By leveraging `SuperSurv`’s screening capabilities, you can scale your
survival analysis to massive genomic and clinical datasets without
sacrificing computational speed or predictive accuracy.

## 5. Identifying the Selected Features

After the Super Learner finishes its cross-validation and fitting, you
likely want to know *which* variables actually survived the screening
process. In clinical and biological research, these selected features
often represent the most important biomarkers or predictive genes.

`SuperSurv` stores this screening information directly inside the fitted
model object, specifically separated into `event.whichScreen` and
`cens.whichScreen`. We can extract the logical matrix that indicates
which columns were passed to each specific base learner.

``` r
# The event.whichScreen object is a logical matrix.
# Rows correspond to the predictors, and columns correspond to the algorithms in our event library.

# 1. Get the logical vector for the 2nd model (Row 2)
is_selected <- fit_highdim$event.whichScreen[2, ] 

# 2. Map it to our original column names
selected_features <- colnames(X_tr)[is_selected]

cat("Total features evaluated:", ncol(X_tr), "\n")
#> Total features evaluated: 29
cat("Features retained by Elastic Net:", length(selected_features), "\n\n")
#> Features retained by Elastic Net: 10

cat("Selected Features:\n")
#> Selected Features:
print(selected_features)
#>  [1] "x0"       "x1"       "x2"       "x3"       "x4"       "x5"      
#>  [7] "x6"       "x7"       "x8"       "noise_10"
```

By extracting these features, you can seamlessly transition from pure
predictive modeling into downstream biological inference or pathway
analysis.
