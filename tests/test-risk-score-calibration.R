library(SuperSurv)

fit_calibration <- getFromNamespace(
  ".fit_risk_score_calibration", "SuperSurv"
)
predict_calibration <- getFromNamespace(
  ".predict_risk_score_survival", "SuperSurv"
)

time <- c(1, 1, 2, 3, 3, 4)
event <- c(1, 1, 1, 0, 1, 0)
score <- c(-0.3, 0.2, 0.5, -0.2, 0.7, 0)
weights <- c(1, 2, 1, 0.5, 1, 1)
centered_score <- score - weighted.mean(score, weights)

breslow <- fit_calibration(
  time, event, score, weights, ties = "breslow"
)
reference_fit <- survival::coxph(
  survival::Surv(time, event) ~ offset(centered_score),
  weights = weights,
  ties = "breslow"
)
reference_baseline <- survival::basehaz(reference_fit, centered = FALSE)
reference_at_events <- vapply(breslow$event.times, function(event_time) {
  tail(reference_baseline$hazard[reference_baseline$time <= event_time], 1L)
}, numeric(1L))
stopifnot(
  max(abs(breslow$cumulative.hazard - reference_at_events)) < 1e-10,
  all(breslow$hazard.increment >= 0)
)

efron <- fit_calibration(time, event, score, weights, ties = "efron")
reference_efron_fit <- survival::coxph(
  survival::Surv(time, event) ~ offset(centered_score),
  weights = weights,
  ties = "efron"
)
reference_efron <- survival::basehaz(reference_efron_fit, centered = FALSE)
reference_efron_at_events <- vapply(efron$event.times, function(event_time) {
  tail(reference_efron$hazard[reference_efron$time <= event_time], 1L)
}, numeric(1L))
stopifnot(
  max(abs(efron$cumulative.hazard - reference_efron_at_events)) < 1e-10
)

evaluation_times <- c(0, 1, 2.5, 5)
prediction <- predict_calibration(
  breslow, score[1:2], evaluation_times, "exponential"
)
shifted <- fit_calibration(
  time, event, score + 100, weights, ties = "breslow"
)
shifted_prediction <- predict_calibration(
  shifted, score[1:2] + 100, evaluation_times, "exponential"
)
product_limit <- predict_calibration(
  breslow, score[1:2], evaluation_times, "product_limit"
)
stopifnot(
  all(prediction[, 1] == 1),
  max(abs(prediction - shifted_prediction)) < 1e-10,
  all(product_limit >= 0 & product_limit <= 1),
  all(apply(product_limit, 1, function(x) all(diff(x) <= 1e-12)))
)

# Validate the stated weighted/unweighted, tied/untied reference matrix.
for (tied in c(FALSE, TRUE)) {
  for (weighted in c(FALSE, TRUE)) {
    for (tie_method in c("breslow", "efron")) {
      reference_time <- if (tied) time else c(1, 1.5, 2, 3, 3.5, 4)
      reference_weights <- if (weighted) weights else rep(1, length(time))
      calibration <- fit_calibration(reference_time, event, score,
                                     reference_weights, ties = tie_method)
      offset_score <- score - calibration$center
      reference <- survival::coxph(
        survival::Surv(reference_time, event) ~ offset_score,
        weights = reference_weights, ties = tie_method,
        init = 1, control = survival::coxph.control(iter.max = 0)
      )
      new_score <- c(-0.2, 0.6)
      reference_survival <- survival::survfit(
        reference, newdata = data.frame(offset_score = new_score - calibration$center),
        stype = 2, ctype = if (tie_method == "efron") 2 else 1
      )
      expected <- t(summary(reference_survival, times = evaluation_times,
                            extend = TRUE)$surv)
      actual <- predict_calibration(calibration, new_score, evaluation_times)
      stopifnot(identical(dim(actual), dim(expected)),
                max(abs(actual - expected)) < 1e-10)
    }
  }
}

# With constant relative risk, Breslow hazard products recover weighted KM.
for (reference_weights in list(rep(1, length(time)), weights)) {
  calibration <- fit_calibration(time, event, rep(0, length(time)),
                                 reference_weights, ties = "breslow")
  actual <- predict_calibration(calibration, 0, evaluation_times, "product_limit")
  reference <- survival::survfit(survival::Surv(time, event) ~ 1,
                                 weights = reference_weights)
  expected <- summary(reference, times = evaluation_times, extend = TRUE)$surv
  stopifnot(max(abs(as.numeric(actual) - expected)) < 1e-10)
}

# Nonconstant relative risks, including factors reaching zero.
hand_calibration <- list(event.times = c(1, 2), hazard.increment = c(0.2, 0.3),
                         cumulative.hazard = c(0.2, 0.5), center = 0)
actual <- predict_calibration(hand_calibration, log(c(1, 2, 6)),
                               c(0, 1, 2, 3), "product_limit")
expected <- rbind(c(1, 0.8, 0.56, 0.56), c(1, 0.6, 0.24, 0.24), c(1, 0, 0, 0))
stopifnot(max(abs(actual - expected)) < 1e-12)

if (requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(20260828)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  glmnet_time <- sample(1:12, n, replace = TRUE)
  glmnet_event <- rbinom(n, 1, 0.65)
  if (sum(glmnet_event) < 10) glmnet_event[seq_len(10)] <- 1

  fit <- surv.glmnet(
    time = glmnet_time,
    event = glmnet_event,
    X = X,
    newdata = X[1:5, , drop = FALSE],
    new.times = c(1, 3, 6, 9),
    obsWeights = rep(1, n),
    id = NULL,
    nfolds = 3,
    ties = "breslow"
  )
  later <- predict(
    fit$fit,
    newdata = X[1:5, , drop = FALSE],
    new.times = c(0, 1, 3, 6, 9, 12)
  )
  stopifnot(
    !is.null(fit$fit$calibration),
    identical(fit$fit$times, fit$fit$calibration$event.times),
    max(abs(fit$pred - later[, c(2, 3, 4, 5), drop = FALSE])) < 1e-10,
    all(later[, 1] == 1),
    all(apply(later, 1, function(x) all(diff(x) <= 1e-12)))
  )
}
