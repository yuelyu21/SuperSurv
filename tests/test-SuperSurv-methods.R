library(SuperSurv)

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]

fit <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X[1:5, , drop = FALSE],
  new.times = seq(20, 100, by = 20),
  event.library = "surv.km",
  cens.library = "surv.km",
  nFolds = 2,
  verbose = FALSE
)

event_w <- event_weights(fit)
censor_w <- censor_weights(fit)
all_weights <- coef(fit, type = "both")
event_learners <- learner_names(fit)
all_learners <- learner_names(fit, type = "both")
times <- eval_times(fit)
model_vars <- training_variables(fit)
selected_vars <- selected_variables(fit, learner = 1)
fit_summary <- summary(fit)

stopifnot(
  inherits(fit, "SuperSurv"),
  is.numeric(event_w),
  is.numeric(censor_w),
  length(event_w) == 1L,
  is.list(all_weights),
  identical(all_weights$event, event_w),
  is.numeric(all_weights$censoring),
  identical(event_learners, names(event_w)),
  identical(all_learners$event, event_learners),
  identical(times, seq(20, 100, by = 20)),
  identical(model_vars, colnames(X)),
  identical(selected_vars, colnames(X)),
  inherits(fit_summary, "summary.SuperSurv"),
  nrow(fit_summary$event) == 1L,
  nrow(fit_summary$censoring) == 1L
)

print(fit)
print(fit_summary)

safe_detect_cores <- getFromNamespace(".safe_detect_cores", "SuperSurv")

stopifnot(
  identical(safe_detect_cores(function() 2), 2L),
  identical(safe_detect_cores(function() NA_integer_), 1L),
  identical(safe_detect_cores(function() NaN), 1L),
  identical(safe_detect_cores(function() Inf), 1L),
  identical(safe_detect_cores(function() 0), 1L),
  identical(safe_detect_cores(function() c(1, 2)), 1L),
  identical(safe_detect_cores(function() "2"), 1L),
  identical(safe_detect_cores(function() stop("boom")), 1L)
)
