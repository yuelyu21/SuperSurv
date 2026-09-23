library(SuperSurv)

expect_error_contains <- function(expr, pattern) {
  message <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
prediction_times <- c(20, 40, 60)

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X[1:5, , drop = FALSE],
  new.times = prediction_times,
  event.library = "surv.km",
  cens.library = "surv.km",
  nFolds = 2,
  verbose = FALSE
)

# Main fitting boundary
expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = as.matrix(X),
    new.times = prediction_times,
    event.library = "surv.km",
    cens.library = "surv.km",
    nFolds = 2
  ),
  "`X` must be a data frame"
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = rep(1, nrow(dat)),
    X = X,
    new.times = prediction_times,
    event.library = "surv.km",
    cens.library = "surv.km",
    nFolds = 2
  ),
  "at least one censored observation"
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    new.times = prediction_times,
    event.library = "surv.km",
    cens.library = "surv.km",
    nFolds = 1
  ),
  "`nFolds` must be one integer between 2"
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    new.times = prediction_times,
    event.library = "surv.km",
    cens.library = "surv.km",
    parallel = c(TRUE, FALSE),
    nFolds = 2
  ),
  "`parallel` must be TRUE or FALSE"
)

# Prediction and object boundaries
expect_error_contains(event_weights(1), "requires a fitted 'SuperSurv' object")
expect_error_contains(learner_names(list()), "requires a fitted 'SuperSurv' object")
expect_error_contains(
  predict(structure(list(), class = "SuperSurv"),
          newdata = X, new.times = prediction_times),
  "missing component(s)"
)
expect_error_contains(
  predict(fit, newdata = as.matrix(X[1:3, ]), new.times = prediction_times),
  "`newdata` must be a data frame"
)
expect_error_contains(
  predict(fit, newdata = X[1:3, -1, drop = FALSE],
          new.times = prediction_times),
  "missing training feature(s)"
)
unexpected <- X[1:3, , drop = FALSE]
unexpected$extra <- 1
expect_error_contains(
  predict(fit, newdata = unexpected, new.times = prediction_times),
  "contains unexpected feature(s)"
)
expect_error_contains(
  predict(fit, newdata = X[1:3, , drop = FALSE], new.times = c(40, 20)),
  "strictly increasing"
)
expect_error_contains(
  predict(fit, newdata = X[1:3, , drop = FALSE],
          new.times = prediction_times, onlySL = 1),
  "`onlySL` must be TRUE or FALSE"
)
expect_error_contains(
  predict(fit, newdata = X[1:3, , drop = FALSE],
          new.times = prediction_times, threshold = -1),
  "`threshold` must be one finite, non-negative number"
)

# Evaluation boundaries
S <- matrix(0.8, nrow = 5, ncol = length(prediction_times))
expect_error_contains(
  eval_brier(dat$duration[1:5], dat$event[1:5], S[-1, ], prediction_times),
  "`S_mat` has dimensions"
)
expect_error_contains(
  eval_logloss(dat$duration[1:5], dat$event[1:5],
               matrix(2, nrow = 5, ncol = 3), prediction_times),
  "survival probabilities in [0, 1]"
)
expect_error_contains(
  eval_cindex(dat$duration[1:5], dat$event[1:5], S, prediction_times,
              eval_time = 40, method = "unknown"),
  "'arg' should be one of"
)
expect_error_contains(
  eval_benchmark(1, X[1:5, , drop = FALSE], dat$duration[1:5],
                 dat$event[1:5], prediction_times),
  "`object` must be a fitted 'SuperSurv' object"
)

# RMST and explanation-plot boundaries
expect_error_contains(
  estimate_marginal_rmst(1, dat, "x0", prediction_times, tau = 40),
  "`fit` must be a fitted 'SuperSurv' object"
)
bad_treatment <- dat
bad_treatment$x0 <- 2
expect_error_contains(
  estimate_marginal_rmst(fit, bad_treatment, "x0", prediction_times, tau = 40),
  "must contain only 0 and 1"
)
expect_error_contains(
  plot_global_importance("not SHAP values"),
  "`shap_values` must be a numeric matrix or data frame"
)
shap_values <- data.frame(x0 = c(0.1, -0.1), x1 = c(0.2, -0.2))
expect_error_contains(
  plot_patient_waterfall(shap_values, patient_index = 3),
  "`patient_index` must be one integer between 1"
)
expect_error_contains(
  plot_dependence(shap_values, data.frame(x0 = 1:2, x1 = 3:4),
                  feature_name = "missing"),
  "present in both `shap_values` and `data`"
)

# Extension contract boundaries
bad_shape_learner <- function(time, event, X, newdata, new.times,
                              obsWeights = NULL, id = NULL, ...) {
  list(pred = matrix(0.5, nrow = 1, ncol = 1), fit = structure(list(), class = "bad_shape"))
}

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    new.times = prediction_times,
    event.library = "bad_shape_learner",
    cens.library = "surv.km",
    nFolds = 2
  ),
  "Learner 'bad_shape_learner' returned an invalid result"
)

bad_screener <- function(X, ...) rep(1, ncol(X))

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    new.times = prediction_times,
    event.library = list(c("surv.km", "bad_screener")),
    cens.library = "surv.km",
    nFolds = 2
  ),
  "Screener 'bad_screener' returned an invalid result"
)
