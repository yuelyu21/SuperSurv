library(SuperSurv)

set.seed(912)
n <- 64L
X <- data.frame(x1 = rnorm(n), x2 = runif(n),
                group = factor(rep(c("A", "B"), length.out = n)))
time <- rexp(n, exp(0.2 * X$x1))
event <- rep(c(1, 1, 1, 0), length.out = n)
newdata <- X[1:4, , drop = FALSE]
times <- c(0, 0.2, 0.6)
expect_error <- function(expr, text) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(text, message, fixed = TRUE))
}
expect_error(surv.coxtime(time, event, X, newdata, times, epochs = 0), "epochs")
expect_error(surv.coxtime(time, event, X, newdata, times, obsWeights = seq_len(n)), "nonuniform")
expect_error(surv.coxtime(time, event, X, newdata, times, reverse = TRUE), "manages")
expect_error(surv.coxtime(time, event, X, newdata, times, standardize_time = NA), "standardize_time")

if (!requireNamespace("survivalmodels", quietly = TRUE)) {
  expect_error(surv.coxtime(time, event, X, newdata, times), "optional R package 'survivalmodels'")
} else if (identical(Sys.getenv("SUPERSURV_RUN_PYCOX_TESTS"), "true")) {
  for (learner in c("surv.deepsurv", "surv.deephit", "surv.coxtime")) {
    args <- list(time = time, event = event, X = X, newdata = newdata,
                 new.times = times, num_nodes = c(8L, 8L), epochs = 2L,
                 batch_size = 16L, batch_norm = FALSE, device = "cpu", seed = 912L)
    if (learner == "surv.deephit") args$cuts <- 8L
    fit <- do.call(learner, args)
    p <- fit$pred
    stopifnot(identical(dim(p), c(4L, 3L)), all(is.finite(p)),
              all(p >= 0 & p <= 1),
              all(apply(p, 1L, function(x) all(diff(x) <= 1e-8))))
    repeated <- predict(fit$fit, newdata[, c("group", "x1", "x2")], new.times = times)
    stopifnot(isTRUE(all.equal(p, repeated, tolerance = 1e-8)))
    expanded <- c(0, 0.1, 0.2, 0.4, 0.6)
    repeated <- predict(fit$fit, newdata, new.times = expanded)
    stopifnot(isTRUE(all.equal(p, repeated[, match(times, expanded)], tolerance = 1e-8)),
              identical(dim(predict(fit$fit, newdata, new.times = 0.2)), c(4L, 1L)),
              identical(dim(predict(fit$fit, newdata[1, , drop = FALSE], new.times = times)), c(1L, 3L)))
    # Compare at the backend's exact time points, including original-scale Cox-Time times.
    model_matrix <- SuperSurv:::.predict_wrapper_matrix(newdata, fit$fit$matrix_spec, learner)
    native <- predict(fit$fit$object, as.data.frame(model_matrix), type = "survival")
    native_times <- as.numeric(colnames(native))
    at_native <- predict(fit$fit, newdata, new.times = native_times)
    stopifnot(max(abs(at_native - native)) < 1e-8)
    expect_error(predict(fit$fit, newdata[, -1], new.times = times), "exactly the features")
    unknown <- newdata
    unknown$group <- factor(rep("unseen", 4))
    expect_error(predict(fit$fit, unknown, new.times = times), "unseen level")
    path <- tempfile(fileext = ".rds")
    saveRDS(fit$fit, path)
    restored <- readRDS(path)
    unlink(path)
    expect_error(predict(restored, newdata, new.times = times), "live Python fitted model")
    args$obsWeights <- seq_len(n)
    expect_error(do.call(learner, args), "nonuniform")
    cat(learner, "end-to-end tests passed.\n")
  }
} else {
  cat("Python fits skipped; set SUPERSURV_RUN_PYCOX_TESTS=true in a configured environment.\n")
}
