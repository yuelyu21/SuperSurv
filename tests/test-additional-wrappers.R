library(SuperSurv)

set.seed(6384)
n <- 60L
X <- data.frame(
  x1 = stats::rnorm(n),
  x2 = stats::runif(n),
  group = factor(rep(c("A", "B"), length.out = n))
)
event_time <- stats::rexp(n, exp(0.4 * X$x1 - 0.3 * X$x2))
censor_time <- stats::rexp(n, 0.35)
time <- pmin(event_time, censor_time)
event <- as.numeric(event_time <= censor_time)
newdata <- X[1:6, , drop = FALSE]
times <- c(0.1, 0.4, 0.8)
weights <- seq(0.5, 1.5, length.out = n)

check_wrapper <- function(result, wrapper_class) {
  stopifnot(
    is.list(result),
    inherits(result$fit, wrapper_class),
    identical(dim(result$pred), c(nrow(newdata), length(times))),
    all(is.finite(result$pred)),
    all(result$pred >= 0 & result$pred <= 1),
    all(apply(result$pred, 1L, function(value) all(diff(value) <= 1e-12)))
  )
  later <- stats::predict(result$fit, newdata = newdata, new.times = times)
  stopifnot(
    identical(dim(later), dim(result$pred)),
    isTRUE(all.equal(later, result$pred, tolerance = 1e-8))
  )
}

if (requireNamespace("mboost", quietly = TRUE)) {
  mboost_fit <- surv.mboost(
    time, event, X, newdata, times, weights, mstop = 20, nu = 0.1
  )
  check_wrapper(mboost_fit, "surv.mboost")
}

if (requireNamespace("flexsurv", quietly = TRUE)) {
  flexsurv_fit <- surv.flexsurvspline(
    time, event, X, newdata, times, weights, k = 1
  )
  check_wrapper(flexsurv_fit, "surv.flexsurvspline")
}

if (requireNamespace("grf", quietly = TRUE)) {
  grf_fit <- surv.grf(
    time, event, X, newdata, times, weights,
    num.trees = 50, min.node.size = 5, seed = 6384
  )
  check_wrapper(grf_fit, "surv.grf")
}

deep_grid <- create_grid(
  "surv.deepsurv", list(epochs = 1L, batch_size = 16L, seed = 6384L)
)
stopifnot(
  inherits(deep_grid, "SuperSurv_grid"),
  length(deep_grid) == 1L,
  exists(deep_grid, mode = "function", inherits = TRUE)
)

capture_error <- function(expression) {
  tryCatch(
    {
      force(expression)
      NA_character_
    },
    error = conditionMessage
  )
}

if (!requireNamespace("survivalmodels", quietly = TRUE)) {
  for (learner in c("surv.deepsurv", "surv.deephit")) {
    arguments <- list(
      time = time, event = event, X = X, newdata = newdata,
      new.times = times, epochs = 1L, batch_size = 16L
    )
    message <- capture_error(do.call(learner, arguments))
    stopifnot(grepl("requires the optional R package 'survivalmodels'", message,
                   fixed = TRUE))
  }
} else if (identical(Sys.getenv("SUPERSURV_RUN_PYCOX_TESTS"), "true")) {
  for (learner in c("surv.deepsurv", "surv.deephit")) {
    arguments <- list(
      time = time, event = event, X = X, newdata = newdata,
      new.times = times, epochs = 1L, batch_size = 16L
    )
    arguments$obsWeights <- weights
    message <- capture_error(do.call(learner, arguments))
    stopifnot(grepl("does not support nonuniform", message, fixed = TRUE))
  }
}
