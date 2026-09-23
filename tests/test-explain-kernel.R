library(SuperSurv)

local({
  expect_missing_backend <- function(explain) {
    error <- tryCatch(
      explain(NULL, data.frame(x = 1), data.frame(x = 1), eval_time = 100),
      error = identity
    )
    stopifnot(
      inherits(error, "error"),
      identical(conditionMessage(error), paste0(
        "Kernel SHAP explanations require the optional 'kernelshap' package. ",
        "Install it to use `explain_kernel()`."
      )),
      is.null(conditionCall(error))
    )
  }

  # Exercise the missing-backend guard even on machines with kernelshap installed.
  explain_without_backend <- explain_kernel
  isolated <- new.env(parent = environment(explain_without_backend))
  isolated$requireNamespace <- function(package, ...) {
    if (identical(package, "kernelshap")) return(FALSE)
    base::requireNamespace(package, ...)
  }
  environment(explain_without_backend) <- isolated
  expect_missing_backend(explain_without_backend)
  if (!requireNamespace("kernelshap", quietly = TRUE)) {
    expect_missing_backend(explain_kernel)
  }
})

if (requireNamespace("kernelshap", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:50, ]
  x_cols <- grep("^x", names(dat))[1:3]
  X <- dat[, x_cols, drop = FALSE]
  wrapper_fit <- surv.coxph(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X[1:3, , drop = FALSE],
    new.times = c(20, 40, 60),
    obsWeights = rep(1, nrow(dat)),
    id = NULL
  )

  explanation <- explain_kernel(
    model = wrapper_fit,
    X_explain = X[1:3, , drop = FALSE],
    X_background = X[4:13, , drop = FALSE],
    nsim = 4, eval_time = 40,
    verbose = FALSE
  )
  reconstructed <- attr(explanation, "baseline") + rowSums(explanation)
  target <- as.numeric(1 - predict(wrapper_fit$fit, X[1:3, ], new.times = 40))
  stopifnot(
    inherits(explanation, "explain"),
    is.data.frame(explanation),
    identical(dim(explanation), c(3L, 3L)),
    identical(names(explanation), names(X)),
    identical(attr(explanation, "backend"), "kernelshap"),
    all(attr(explanation, "converged")),
    max(abs(reconstructed - attr(explanation, "predictions"))) < 1e-6,
    max(abs(reconstructed - target)) < 1e-6,
    identical(attr(explanation, "target"), "event_probability")
  )

  set.seed(6384)
  ensemble <- SuperSurv(
    time = dat$duration, event = dat$event, X = X,
    event.library = c("surv.coxph", "surv.weibull", "surv.km"),
    cens.library = "surv.km", nFolds = 3,
    control = list(saveFitLibrary = TRUE), verbose = FALSE
  )
  # Force a heterogeneous mixture, including a weight below the old cutoff.
  ensemble$event.coef[] <- c(0.45, 0.5495, 0.0005)
  for (horizon in c(20, 40)) {
    mixed <- explain_kernel(ensemble, X[1:3, ], X[4:13, ],
                            nsim = 4, eval_time = horizon)
    target <- as.numeric(1 - predict(ensemble, X[1:3, ], horizon, type = "event"))
    background <- as.numeric(1 - predict(ensemble, X[4:13, ], horizon, type = "event"))
    stopifnot(
      max(abs(attr(mixed, "predictions") - target)) < 1e-10,
      abs(attr(mixed, "baseline") - mean(background)) < 1e-10,
      max(abs(attr(mixed, "baseline") + rowSums(mixed) - target)) < 1e-6
    )
  }
  best <- explain_kernel(ensemble, X[1:3, ], X[4:13, ],
                         only_best = TRUE, eval_time = 40)
  best_target <- as.numeric(1 - predict(ensemble$event.fitLibrary[[2]],
                                       X[1:3, ], new.times = 40))
  stopifnot(max(abs(attr(best, "predictions") - best_target)) < 1e-10)
  for (bad_time in list(NULL, NA_real_, Inf, -1, c(1, 2))) {
    error <- tryCatch(explain_kernel(ensemble, X[1:3, ], X[4:13, ],
                                     eval_time = bad_time),
                      error = conditionMessage)
    stopifnot(is.character(error), grepl("eval_time", error, fixed = TRUE))
  }
}
