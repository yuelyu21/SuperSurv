library(SuperSurv)

set.seed(9126384)
n <- 120L
X <- data.frame(x1 = rnorm(n), x2 = runif(n),
                group = factor(rep(c("A", "B"), length.out = n)))
failure <- rexp(n, exp(0.3 * X$x1 - 0.2 * X$x2))
censor <- rexp(n, 0.25)
time <- pmin(failure, censor)
event <- as.numeric(failure <= censor)
newdata <- X[1:5, , drop = FALSE]
times <- c(0, 0.2, 0.6, 1.2)
weights <- seq(0.7, 1.3, length.out = n)

expect_error <- function(expr, text) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(text, message, fixed = TRUE))
}

check_fit <- function(result) {
  p <- result$pred
  stopifnot(identical(dim(p), c(5L, 4L)), all(is.finite(p)),
            all(p >= 0 & p <= 1), all(p[, 1] == 1),
            all(apply(p, 1L, function(x) all(diff(x) <= 1e-8))))
  later <- predict(result$fit, newdata[, c("group", "x2", "x1")], new.times = times)
  stopifnot(isTRUE(all.equal(p, later, tolerance = 1e-10)))
  expanded <- sort(unique(c(times, 0.4, 2)))
  later <- predict(result$fit, newdata, new.times = expanded)
  stopifnot(isTRUE(all.equal(p, later[, match(times, expanded)], tolerance = 1e-10)))
  one_time <- predict(result$fit, newdata, new.times = 0.6)
  one_person <- predict(result$fit, newdata[1, , drop = FALSE], new.times = times)
  one_cell <- predict(result$fit, newdata[1, , drop = FALSE], new.times = 0.6)
  stopifnot(identical(dim(one_time), c(5L, 1L)),
            identical(dim(one_person), c(1L, 4L)),
            identical(dim(one_cell), c(1L, 1L)),
            max(abs(one_time[, 1] - p[, 3])) < 1e-10,
            max(abs(one_person[1, ] - p[1, ])) < 1e-10)
  path <- tempfile(fileext = ".rds")
  saveRDS(result$fit, path)
  restored <- readRDS(path)
  unlink(path)
  stopifnot(isTRUE(all.equal(p, predict(restored, newdata, new.times = times),
                            tolerance = 1e-10)))
  unknown <- newdata
  unknown$group <- factor(rep("unseen", nrow(unknown)))
  expect_error(predict(result$fit, unknown, new.times = times), "unseen level")
  expect_error(predict(result$fit, newdata[, -1], new.times = times), "exactly the features")
  expect_error(predict(result$fit, newdata, new.times = c(1, 0)), "new.times")
}

expect_error(surv.flexsurvreg(time, event, X, newdata, times, dist = "bad"), "`dist`")
expect_error(surv.flexsurvreg(time, event, X, newdata, times, weights = weights), "manages")
expect_error(surv.survPen(time, event, X, newdata, times, obsWeights = weights), "nonuniform")
expect_error(surv.survPen(time, event, X, newdata, times, baseline.df = 1), "baseline.df")
expect_error(surv.survPen(time, event, X, newdata, times, n.legendre = NA), "n.legendre")
expect_error(surv.survPen(time, event, X, newdata, times, expected = rep(0, n)), "manages")
expect_error(surv.survPen(time, event, X, newdata, times, formula = time ~ x1), "one-sided")

if (requireNamespace("flexsurv", quietly = TRUE)) {
  for (distribution in c("gengamma", "gompertz", "gamma", "weibull", "exp", "lnorm", "llogis")) {
    result <- surv.flexsurvreg(time, event, X, newdata, times, weights,
                               dist = distribution)
    check_fit(result)
    native <- summary(result$fit$object, newdata = newdata, t = times,
                       type = "survival", ci = FALSE)
    native <- do.call(rbind, lapply(native, function(x) x$est))
    stopifnot(max(abs(result$pred - native)) < 1e-8)
  }
  # Regression for the shared flexsurv singleton-time extraction helper.
  spline <- surv.flexsurvspline(time, event, X, newdata, times, weights, k = 0)
  stopifnot(identical(dim(predict(spline$fit, newdata, new.times = 0.6)), c(5L, 1L)))
} else {
  expect_error(surv.flexsurvreg(time, event, X, newdata, times), "optional package 'flexsurv'")
}

if (requireNamespace("survPen", quietly = TRUE)) {
  fits <- list(
    surv.survPen(time, event, X, newdata, times, baseline.df = 3),
    surv.survPen(time, event, X, newdata, times,
                 formula = ~ smf(.supersurv_time, df = 3) + smf(x1, df = 3) + x2 + group),
    surv.survPen(time, event, X, newdata, times,
                 formula = ~ tensor(.supersurv_time, x1, df = c(3, 3)) + x2 + group)
  )
  for (result in fits) {
    check_fit(result)
    for (j in seq_along(times)) {
      native_data <- newdata
      native_data$.supersurv_time <- times[j]
      native <- predict(result$fit$object, native_data, n.legendre = 50L)$surv
      stopifnot(max(abs(result$pred[, j] - native)) < 1e-8)
    }
  }
  expect_error(surv.survPen(time, event, X, newdata, times, formula = ~ absent),
               "not in the training features")
  expect_error(predict(fits[[1]]$fit, newdata, new.times = times, n.legendre = 0),
               "n.legendre")
} else {
  expect_error(surv.survPen(time, event, X, newdata, times), "optional package 'survPen'")
}

cat("Flexible-wrapper validation, native agreement, prediction, and persistence tests passed.\n")
