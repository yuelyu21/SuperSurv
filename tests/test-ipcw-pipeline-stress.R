library(SuperSurv)

# Deliberately misspecified, fixed survival curves force active truncation.
# These fixtures test pipeline mechanics, not statistical threshold selection.
stress_state <- new.env(parent = emptyenv())
stress_state$calls <- list()
stress_state$censor_scale <- 1
make_stress_learner <- function(name, rate, censoring = FALSE) {
  force(name)
  force(rate)
  force(censoring)
  function(time, event, X, newdata, new.times, ...) {
    stress_state$calls[[length(stress_state$calls) + 1L]] <- list(
      learner = name, train = X$row_id, predict = newdata$row_id
    )
    actual_rate <- rate * if (censoring) stress_state$censor_scale else 1
    pred <- matrix(rep(exp(-actual_rate * new.times), each = nrow(newdata)),
                   nrow = nrow(newdata))
    list(pred = pred, fit = list(rate = actual_rate))
  }
}
surv.stress_event_slow <- make_stress_learner("event_slow", 0.2)
surv.stress_event_fast <- make_stress_learner("event_fast", 1)
surv.stress_censor_slow <- make_stress_learner("censor_slow", 8, TRUE)
surv.stress_censor_fast <- make_stress_learner("censor_fast", 10, TRUE)
surv.stress_initializer <- make_stress_learner("initializer", 8, TRUE)

stress_n <- 60L
stress_X <- data.frame(row_id = seq_len(stress_n), x = rep(c(-1, 1), 30))
stress_time <- seq(0.1, 3, length.out = stress_n)
stress_event <- rep(c(1, 0), 30)
stress_grid <- c(0, 0.5, 1, 1.5, 2, 2.5)
stress_folds <- split(seq_len(stress_n), rep(seq_len(3), length.out = stress_n))
stress_newX <- data.frame(row_id = 101:104, x = c(-1, 1, -1, 1))

run_stress <- function(method, floor = 1e-4, cap = 100,
                       censor_scale = 1, warn_fraction = 0.05, max_iter = 20) {
  stress_state$calls <- list()
  stress_state$censor_scale <- censor_scale
  warnings <- character()
  fit <- withCallingHandlers(
    SuperSurv(
      time = stress_time, event = stress_event, X = stress_X,
      newdata = stress_newX, new.times = stress_grid,
      event.library = c("surv.stress_event_slow", "surv.stress_event_fast"),
      cens.library = c("surv.stress_censor_slow", "surv.stress_censor_fast"),
      metalearner = method, nFolds = 3,
      cvControl = list(validRows = stress_folds),
      control = list(
        initWeightAlg = "surv.stress_initializer",
        event.t.grid = stress_grid, cens.t.grid = stress_grid,
        ipcw.floor = floor, ipcw.cap = cap, max.SL.iter = max_iter,
        truncation.warn.fraction = warn_fraction, saveFitLibrary = FALSE
      ),
      verbose = FALSE, parallel = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warnings, calls = stress_state$calls)
}

check_routing <- function(result) {
  calls <- result$calls
  # Four candidates x three folds, one full-sample initializer, four refits.
  stopifnot(length(calls) == 17L, calls[[13]]$learner == "initializer")
  for (call in calls[seq_len(12)]) {
    stopifnot(length(call$train) == 40L, length(call$predict) == 20L,
              length(intersect(call$train, call$predict)) == 0L,
              setequal(c(call$train, call$predict), stress_X$row_id))
  }
  stopifnot(identical(calls[[13]]$train, stress_X$row_id),
            identical(calls[[13]]$predict, stress_X$row_id))
  for (call in calls[14:17]) {
    stopifnot(identical(call$train, stress_X$row_id),
              identical(call$predict, stress_newX$row_id))
  }
}

check_fit <- function(result, cap, active = TRUE) {
  fit <- result$fit
  d <- fit$algorithm.diagnostics
  check_routing(result)
  stopifnot(d$converged, d$event.optimizer$converged,
            d$censoring.optimizer$converged,
            is.finite(d$event.optimizer$objective),
            is.finite(d$censoring.optimizer$objective),
            all(is.finite(c(fit$event.cvRisks, fit$cens.cvRisks))),
            !any(c(fit$event.errorsInCVLibrary, fit$cens.errorsInCVLibrary,
                   fit$event.errorsInLibrary, fit$cens.errorsInLibrary)),
            d$ipcw.maximum.stabilized.weight <= cap + 1e-10,
            is.finite(d$ipcw.minimum.effective.sample.size),
            d$ipcw.minimum.effective.sample.size > 0)
  for (weights in list(fit$event.coef, fit$cens.coef)) {
    stopifnot(all(is.finite(weights)), all(weights >= 0),
              abs(sum(weights) - 1) < 1e-8)
  }
  for (pred in list(fit$event.predict, fit$cens.predict)) {
    stopifnot(identical(dim(pred), c(4L, length(stress_grid))),
              # Convex sums can overshoot a boundary by floating-point roundoff.
              all(is.finite(pred)), all(pred >= -1e-12 & pred <= 1 + 1e-12),
              all(apply(pred, 1, function(x) all(diff(x) <= 1e-10))))
  }
  if (active) {
    stopifnot(d$ipcw.minimum.denominator < 1e-4,
              d$ipcw.maximum.raw.weight > cap,
              abs(d$ipcw.maximum.stabilized.weight - cap) < 1e-8,
              d$ipcw.floor.fraction > 0, d$ipcw.cap.fraction > 0,
              sum(grepl("IPCW stabilization affected", result$warnings)) == 1L)
  } else {
    stopifnot(d$ipcw.floor.fraction == 0, d$ipcw.cap.fraction == 0,
              !any(grepl("IPCW stabilization affected", result$warnings)))
  }
  stopifnot(!any(grepl("did not converge|failed", result$warnings)))
}

check_event_reference <- function(result) {
  fit <- result$fit
  floor <- fit$control$ipcw.floor
  cap <- fit$control$ipcw.cap
  # Reconstruct the final event update without using package IPCW helpers.
  G_grid <- drop(cbind(exp(-8 * stress_grid), exp(-10 * stress_grid)) %*%
                   fit$cens.coef)
  G_left <- vapply(stress_time, function(t) {
    before <- which(stress_grid < t)
    if (length(before)) G_grid[max(before)] else 1
  }, numeric(1))
  time_long <- rep(stress_time, length(stress_grid))
  event_long <- rep(stress_event, length(stress_grid))
  grid_long <- rep(stress_grid, each = stress_n)
  G_left_long <- rep(G_left, length(stress_grid))
  G_grid_long <- rep(G_grid, each = stress_n)
  failure <- time_long <= grid_long & event_long == 1
  survivor <- time_long > grid_long
  inv_G_left <- pmin(1 / pmax(G_left_long, floor), cap)
  inv_G_grid <- pmin(1 / pmax(G_grid_long, floor), cap)
  S_grid <- drop(cbind(exp(-0.2 * stress_grid), exp(-stress_grid)) %*%
                   fit$event.coef)
  S_long <- rep(S_grid, each = stress_n)
  d <- fit$algorithm.diagnostics$event.optimizer
  if (d$method == "logloss") {
    probability <- pmin(pmax(S_long, fit$control$logloss.eps),
                        1 - fit$control$logloss.eps)
    contribution <- -failure * inv_G_left * log(1 - probability) -
      survivor * inv_G_grid * log(probability)
    active <- failure | survivor
    denominator <- ifelse(failure, G_left_long, G_grid_long)[active]
    diagnostic <- d$combined.weights
  } else {
    contribution <- (1 - failure * inv_G_left - S_long)^2
    active <- failure
    denominator <- G_left_long[active]
    diagnostic <- d$weighted.rows
  }
  inverse <- 1 / pmax(denominator, floor)
  weight <- pmin(inverse, cap)
  stopifnot(
    isTRUE(all.equal(d$objective, mean(contribution), tolerance = 1e-10)),
    diagnostic$rows == sum(active),
    diagnostic$floor.count == sum(denominator < floor),
    diagnostic$cap.count == sum(inverse > cap),
    isTRUE(all.equal(diagnostic$effective.sample.size,
                     sum(weight)^2 / sum(weight^2), tolerance = 1e-12))
  )
  if (d$method == "logloss") {
    stopifnot(isTRUE(all.equal(
      diagnostic$capped.loss.fraction,
      sum(contribution[active][inverse > cap]) / sum(contribution[active]),
      tolerance = 1e-12
    )))
  }
}

stress_results <- list()
for (method in c("brier", "logloss")) {
  for (cap in c(200, 100, 50)) {
    result <- run_stress(method, cap = cap)
    check_fit(result, cap)
    check_event_reference(result)
    stress_results[[paste(method, cap, sep = "_")]] <- result
  }
  # Equal effective thresholds via different floor/cap settings must agree.
  floor_result <- run_stress(method, floor = 0.02, cap = 100)
  cap_result <- stress_results[[paste(method, 50, sep = "_")]]
  check_routing(floor_result)
  check_event_reference(floor_result)
  stopifnot(floor_result$fit$algorithm.diagnostics$converged,
            isTRUE(all.equal(floor_result$fit$event.coef,
                             cap_result$fit$event.coef, tolerance = 1e-10)),
            isTRUE(all.equal(floor_result$fit$cens.coef,
                             cap_result$fit$cens.coef, tolerance = 1e-10)),
            isTRUE(all.equal(floor_result$fit$event.predict,
                             cap_result$fit$event.predict, tolerance = 1e-10)))
  quiet <- run_stress(method, warn_fraction = 1)
  stopifnot(!any(grepl("IPCW stabilization affected", quiet$warnings)),
            quiet$fit$algorithm.diagnostics$ipcw.cap.fraction > 0)
  inactive <- run_stress(method, censor_scale = 0.01)
  check_fit(inactive, cap = 100, active = FALSE)
  limited <- run_stress(method, max_iter = 1)
  stopifnot(!limited$fit$algorithm.diagnostics$converged,
            any(grepl("did not converge", limited$warnings)),
            any(grepl("IPCW stabilization affected", limited$warnings)))
}

stress_summary <- do.call(rbind, lapply(names(stress_results), function(name) {
  result <- stress_results[[name]]
  d <- result$fit$algorithm.diagnostics
  data.frame(setting = name, iterations = d$iterations,
             threshold = d$ipcw.effective.denominator.threshold,
             minimum_denominator = d$ipcw.minimum.denominator,
             max_raw_weight = d$ipcw.maximum.raw.weight,
             max_stabilized_weight = d$ipcw.maximum.stabilized.weight,
             floor_fraction = d$ipcw.floor.fraction,
             cap_fraction = d$ipcw.cap.fraction,
             minimum_ess = d$ipcw.minimum.effective.sample.size,
             converged = d$converged,
             event_objective = d$event.optimizer$objective)
}))
print(stress_summary, row.names = FALSE)
cat("PASS: active/inactive truncation, warning controls, threshold equivalence,",
    "nonconvergence reporting, and fold/refit routing.\n")
