# 6. Parametric Survival Models

## Introduction

While the Cox model is semi-parametric (it leaves the baseline hazard
unspecified), fully parametric models assume that survival times follow
a specific mathematical distribution, such as the Weibull, Exponential,
or Log-Normal distribution.

Parametric models are statistically powerful because they produce
perfectly smooth survival curves. However, they are highly brittle: if
you assume the data follows a Weibull distribution, but the true
biological hazard has a completely different shape, the model will be
heavily biased.

`SuperSurv` acts as a **safety net**. You can include multiple
parametric assumptions in your library. If a parametric assumption
perfectly matches your data, `SuperSurv` will give it a high weight. If
the assumption is wrong, the cross-validation risk will spike, and
`SuperSurv` will safely assign it a weight of zero.

## 1. Setup and Library Definition

``` r
library(SuperSurv)
library(survival)

data("metabric", package = "SuperSurv")
set.seed(42)

train_idx <- sample(1:nrow(metabric), 0.7 * nrow(metabric))
train <- metabric[train_idx, ]
test  <- metabric[-train_idx, ]

X_tr <- train[, grep("^x", names(metabric))]
X_te <- test[, grep("^x", names(metabric))]
new.times <- seq(50, 200, by = 25)

# Define a library covering different parametric assumptions
parametric_library <- c("surv.coxph",       # Semi-parametric baseline
                        "surv.weibull",     # Assumes hazard increases/decreases monotonically
                        "surv.exponential", # Assumes constant hazard over time
                        "surv.lognormal")   # Assumes hazard rises then falls
```

## 2. Fitting the Parametric Ensemble

We run the ensemble exactly as before. Internally, `SuperSurv` will fit
these Accelerated Failure Time (AFT) models and map their continuous
survival predictions onto our discrete `new.times` evaluation grid.

``` r
fit_parametric <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newX = X_te,
  new.times = new.times,
  event.SL.library = parametric_library,
  cens.SL.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE),
  verbose = FALSE,
  selection = "ensemble",
  nFolds = 3
)
```

## 3. Evaluating the “Safety Net”

Let’s look at the cross-validated risks and the final meta-learner
weights.

``` r
cat("Cross-Validated Risks (Lower is better):\n")
#> Cross-Validated Risks (Lower is better):
print(round(fit_parametric$event.cvRisks, 4))
#>       surv.coxph_screen.all     surv.weibull_screen.all 
#>                      8.8384                      8.8420 
#> surv.exponential_screen.all   surv.lognormal_screen.all 
#>                      8.8508                      8.8444

cat("\nFinal Ensemble Weights:\n")
#> 
#> Final Ensemble Weights:
print(round(fit_parametric$event.coef, 4))
#>       surv.coxph_screen.all     surv.weibull_screen.all 
#>                      0.9979                      0.0000 
#> surv.exponential_screen.all   surv.lognormal_screen.all 
#>                      0.0000                      0.0021
```

### Interpretation

Look closely at the weights assigned to `surv.exponential`. The
Exponential distribution assumes that the risk of the event (hazard) is
completely constant over time. In real-world cancer datasets like
`metabric`, this assumption is almost always false (risk usually
increases with time or peaks shortly after surgery).

Because the Exponential assumption fits the data poorly, its
cross-validated risk will be high, and `SuperSurv` will smartly assign
it a weight of $0.00$.

By including parametric models in your `SuperSurv` library, you allow
the data—not the researcher—to dictate which mathematical distributions
are actually appropriate for your patient cohort!
