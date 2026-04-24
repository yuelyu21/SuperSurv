## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
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

## ----define-library-----------------------------------------------------------
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

## ----fit-screen, results='hide', message=FALSE, warning=FALSE-----------------
fit_highdim <- SuperSurv(
  time = train$duration,
  event = train$event,
  X = X_tr,
  newdata = X_te,
  new.times = new.times,
  event.library = my_screen_library, # Pass our screened library
  cens.library = cens_library,
  control = list(saveFitLibrary = TRUE),
  verbose = T,
  nFolds = 3
)

## ----evaluate-screen----------------------------------------------------------
# Check the integrated performance
screened_performance <- eval_summary(
  object = fit_highdim,
  newdata = X_te,
  time = test$duration,
  event = test$event,
  eval_times = new.times
)

# Plot the benchmark
# plot_benchmark(fit_highdim, newdata = X_te, time = test$duration, event = test$event, eval_times = new.times)

## ----extract-features---------------------------------------------------------
# Get the variables retained by the second event learner.
selected_features <- selected_variables(fit_highdim, learner = 2)

cat("Total features evaluated:", ncol(X_tr), "\n")
cat("Features retained by Elastic Net:", length(selected_features), "\n\n")

cat("Selected Features:\n")
print(selected_features)

