library(SuperSurv)
library(survival)

# Exact event-time ties exercise G(T-) rather than G(T).
time <- c(1, 1, 2, 2, 3, 4)
event <- c(1, 0, 1, 0, 1, 0)
times <- c(1, 2, 3)
survival_prediction <- matrix(
  c(1, 0.8, 0.5, 0, 0.7, 0.4),
  nrow = length(time), ncol = length(times), byrow = TRUE
)
survival_prediction <- survival_prediction[, rep(seq_len(ncol(survival_prediction)),
                                                  each = 2), drop = FALSE]
survival_prediction <- survival_prediction[, c(1, 3, 5), drop = FALSE]

reverse_grid <- sort(unique(time))
reverse_survival <- cumprod(vapply(reverse_grid, function(value) {
  at_risk <- sum(time >= value)
  failures <- sum(time == value & event == 1)
  censorings <- sum(time == value & event == 0)
  if (censorings == 0) 1 else 1 - censorings / (at_risk - failures)
}, numeric(1L)))
step_at <- function(x, left = FALSE) {
  vapply(x, function(value) {
    index <- which(if (left) reverse_grid < value else reverse_grid <= value)
    if (length(index)) reverse_survival[max(index)] else 1
  }, numeric(1L))
}
expected <- vapply(seq_along(times), function(j) {
  horizon <- times[j]
  weight <- ifelse(
    time <= horizon & event == 1, 1 / step_at(time, left = TRUE),
    ifelse(time > horizon, 1 / step_at(horizon), 0)
  )
  mean(weight * (as.integer(time > horizon) - survival_prediction[, j])^2)
}, numeric(1L))
observed <- eval_brier(time, event, survival_prediction, times)$brier_scores
stopifnot(isTRUE(all.equal(observed, expected, tolerance = 1e-12)))

# Brier loss accepts valid probability boundaries without log-loss clipping.
boundary <- eval_brier(
  time = c(1, 2), event = c(1, 1),
  S_mat = matrix(c(0, 0, 1, 0), nrow = 2, byrow = TRUE),
  times = c(1, 2)
)
stopifnot(all(is.finite(boundary$brier_scores)))

boundary_logloss <- SuperSurv:::eval_logloss(
  time = c(1, 2), event = c(1, 1),
  S_mat = matrix(c(0, 0, 1, 0), nrow = 2, byrow = TRUE),
  times = c(1, 2), eps = 1e-8
)
stopifnot(
  all(is.finite(boundary_logloss$logloss_scores)),
  identical(
    boundary_logloss$diagnostics$censoring.source,
    "marginal reverse Kaplan-Meier"
  ),
  boundary_logloss$diagnostics$probability.clip.fraction > 0
)

# The backward-compatible default reproduces the original marginal-KM formula.
expected_logloss <- vapply(seq_along(times), function(j) {
  horizon <- times[j]
  probability <- pmin(pmax(survival_prediction[, j], 1e-10), 1 - 1e-10)
  failure_weight <- as.numeric(time <= horizon & event == 1) /
    step_at(time, left = TRUE)
  survivor_weight <- as.numeric(time > horizon) / step_at(horizon)
  mean(-failure_weight * log(1 - probability) -
         survivor_weight * log(probability))
}, numeric(1L))
observed_logloss <- eval_logloss(
  time, event, survival_prediction, times
)
stopifnot(isTRUE(all.equal(
  observed_logloss$logloss_scores, expected_logloss, tolerance = 1e-12
)))
positional_logloss <- eval_logloss(
  time, event, survival_prediction, times,
  min(times), max(times), 1e-6, 1e-8
)
named_logloss <- eval_logloss(
  time, event, survival_prediction, times,
  tmin = min(times), tmax = max(times), ipcw_floor = 1e-6, eps = 1e-8
)
stopifnot(isTRUE(all.equal(positional_logloss, named_logloss)))

conditional_logloss <- eval_logloss(
  time = c(1, 2, 3, 4),
  event = c(1, 0, 1, 0),
  S_mat = matrix(c(
    0.9, 0.7, 0.5,
    0.8, 0.6, 0.4,
    0.7, 0.5, 0.3,
    0.6, 0.4, 0.2
  ), nrow = 4, byrow = TRUE),
  times = c(1, 2, 3),
  ipcw_floor = 1e-4,
  ipcw_cap = 50,
  G_T_left = c(1, 0.5, 0.01, 0.25),
  G_times = matrix(c(
    1, 0.8, 0.6,
    1, 0.5, 0.01,
    1, 0.4, 0.01,
    1, 0.3, 0.01
  ), nrow = 4, byrow = TRUE)
)
conditional_diagnostics <- conditional_logloss$diagnostics
stopifnot(
  all(is.finite(conditional_logloss$logloss_scores)),
  identical(conditional_diagnostics$censoring.source, "user supplied"),
  identical(conditional_diagnostics$effective.denominator.threshold, 0.02),
  conditional_diagnostics$combined.weights$cap.count > 0,
  conditional_diagnostics$combined.weights$effective.sample.size > 0,
  conditional_diagnostics$combined.weights$capped.loss.fraction >= 0,
  conditional_diagnostics$combined.weights$capped.loss.fraction <= 1
)

missing_G_error <- tryCatch(
  eval_logloss(
    time = c(1, 2), event = c(1, 0),
    S_mat = matrix(c(0.8, 0.6), nrow = 2),
    times = 1, G_T_left = c(1, 0.5)
  ),
  error = conditionMessage
)
invalid_G_error <- tryCatch(
  eval_logloss(
    time = c(1, 2), event = c(1, 0),
    S_mat = matrix(c(0.8, 0.6), nrow = 2),
    times = 1, G_T_left = c(1, 0.5), G_times = matrix(c(1, 1.2), nrow = 2)
  ),
  error = conditionMessage
)
stopifnot(
  grepl("Supply both", missing_G_error, fixed = TRUE),
  grepl("finite probability matrix", invalid_G_error, fixed = TRUE)
)

predict.SuperSurv_benchmark_test <- function(object, newdata, new.times, ...) {
  object$predictions
}
ensemble <- matrix(c(0.9, 0.7, 0.4, 0.8, 0.6, 0.3,
                     0.7, 0.4, 0.2, 0.5, 0.3, 0.1), nrow = 4, byrow = TRUE)
library_array <- array(
  c(ensemble, pmin(ensemble + 0.05, 1)),
  dim = c(4, 3, 2),
  dimnames = list(NULL, NULL, c("learner_a", "learner_b"))
)
mock <- structure(
  list(predictions = list(
    event.predict = ensemble,
    event.library.predict = library_array
  )),
  class = "SuperSurv_benchmark_test"
)
benchmark <- SuperSurv:::eval_benchmark(
  mock, newdata = data.frame(x = 1:4),
  time = c(1, 2, 3, 4), event = c(1, 1, 0, 1), eval_times = c(1, 2, 3)
)
stopifnot(
  inherits(benchmark, "SuperSurv_benchmark"),
  identical(names(benchmark), c("summary", "by_time", "eval_times", "risk_time")),
  nrow(benchmark$summary) == 3L,
  nrow(benchmark$by_time) == 9L,
  "IPCW_LogLoss" %in% names(benchmark$summary),
  "LogLoss" %in% names(benchmark$by_time)
)

summary_result <- eval_summary(
  mock, newdata = data.frame(x = 1:4),
  time = c(1, 2, 3, 4), event = c(1, 1, 0, 1), eval_times = c(1, 2, 3)
)
stopifnot(
  inherits(summary_result, "SuperSurv_eval"),
  is.data.frame(attr(summary_result, "by_time"))
)

fallback <- SuperSurv:::.combine_benchmark_plots(
  list(brier = structure(list(), class = "mock_plot"),
       auc = structure(list(), class = "mock_plot")),
  patchwork_available = FALSE
)
stopifnot(inherits(fallback, "SuperSurv_plot_list"), length(fallback) == 2L)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  # No prediction method is available for this object: plotting must reuse it.
  benchmark_plot <- plot_benchmark(benchmark, metrics = "brier")
  stopifnot(inherits(benchmark_plot, "ggplot"),
            identical(benchmark_plot$data$Brier, benchmark$by_time$Brier))
  malformed <- structure(list(by_time = data.frame()), class = "SuperSurv_benchmark")
  stopifnot(grepl("by_time", tryCatch(
    plot_benchmark(malformed), error = conditionMessage), fixed = TRUE))
  stopifnot(grepl("omit", tryCatch(
    plot_benchmark(benchmark, eval_times = 1), error = conditionMessage), fixed = TRUE))
}

if (requireNamespace("riskRegression", quietly = TRUE)) {
  lung <- survival::lung
  lung_time <- lung$time
  lung_event <- as.integer(lung$status == 2)
  score_times <- c(100, 300, 500)
  km_fit <- survival::survfit(survival::Surv(lung_time, lung_event) ~ 1)
  km_survival <- summary(km_fit, times = score_times, extend = TRUE)$surv
  prediction <- matrix(
    km_survival, nrow = length(lung_time), ncol = length(score_times),
    byrow = TRUE
  )
  supersurv_score <- eval_brier(
    lung_time, lung_event, prediction, score_times
  )$brier_scores
  reference <- riskRegression::Score(
    object = list(KM = 1 - prediction),
    formula = Surv(time, event) ~ 1,
    data = data.frame(time = lung_time, event = lung_event),
    times = score_times, metrics = "brier", cens.model = "km",
    summary = "risk", conf.int = FALSE, se.fit = FALSE
  )
  reference_score <- reference$Brier$score
  reference_score <- reference_score[reference_score$model == "KM", "Brier"][[1L]]
  stopifnot(max(abs(supersurv_score - reference_score)) < 1e-6)

  tied_time <- rep(1:8, each = 8)
  tied_event <- rep(c(1, 0, 1, 0, 1, 1, 0, 0), 8)
  tied_times <- c(2, 4, 6)
  tied_prediction <- outer(
    seq_along(tied_time), tied_times,
    function(i, t) pmax(0.05, pmin(0.95, exp(-t * (0.03 + i / 500))))
  )
  tied_supersurv <- eval_brier(
    tied_time, tied_event, tied_prediction, tied_times
  )$brier_scores
  tied_reference <- riskRegression::Score(
    object = list(model = 1 - tied_prediction),
    formula = Surv(time, event) ~ 1,
    data = data.frame(time = tied_time, event = tied_event),
    times = tied_times, metrics = "brier", cens.model = "km",
    summary = "risk", conf.int = FALSE, se.fit = FALSE
  )$Brier$score
  tied_reference <- tied_reference[tied_reference$model == "model", "Brier"][[1L]]
  stopifnot(max(abs(tied_supersurv - tied_reference)) < 1e-6)
}
