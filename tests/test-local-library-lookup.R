library(SuperSurv)

set.seed(638422)
n <- 72L
X <- data.frame(x1 = rnorm(n), x2 = runif(n))
failure <- rexp(n, exp(0.2 * X$x1))
censoring <- rexp(n, 0.6)
time <- pmin(failure, censoring)
event <- as.numeric(failure <= censoring)
times <- as.numeric(quantile(time, c(0.15, 0.35, 0.55)))
grid <- c(0, times)
folds <- lapply(1:3, function(v) which(rep(1:3, length.out = n) == v))
common <- list(
  time = time, event = event, X = X, newdata = X[1:4, ], new.times = times,
  nFolds = 3L, cvControl = list(validRows = folds),
  control = list(event.t.grid = grid, cens.t.grid = grid, saveFitLibrary = TRUE)
)

local_fit <- function(parallel = FALSE, method = "brier") {
  calls <- new.env(parent = emptyenv())
  calls$learners <- 0L
  calls$screens <- 0L
  surv.local_scope_base <- function(time, event, X, newdata, new.times,
                                   obsWeights = NULL, id = NULL, marker, ...) {
    stopifnot(ncol(X) == ncol(newdata), marker %in% c(17L, 19L))
    # CV also evaluates unused learner/screener combinations. Count the
    # requested one-column combinations, plus their final refits, only.
    if (ncol(X) == 1L) calls$learners <- calls$learners + 1L
    out <- SuperSurv::surv.km(time, event, X, newdata, new.times, obsWeights, id)
    out$fit$local_scope_marker <- marker
    out
  }
  surv.local_scope_direct <- function(...) SuperSurv::surv.coxph(...)
  screen.local_scope <- function(X, ...) {
    calls$screens <- calls$screens + 1L
    seq_len(ncol(X)) == 1L
  }
  learners <- create_grid("surv.local_scope_base", list(marker = c(17L, 19L)))
  # Each grid member retains the selected base function, not a later rebinding.
  surv.local_scope_base <- function(...) stop("The grid used a rebound base function.")
  args <- c(common, list(
    event.library = list(c(learners[1L], "screen.local_scope"), "surv.local_scope_direct"),
    cens.library = list(c(learners[2L], "screen.local_scope")),
    parallel = parallel, metalearner = method
  ))
  fit <- do.call(SuperSurv, args)
  if (!parallel) stopifnot(calls$learners == 8L, calls$screens == 8L)
  stopifnot(
    inherits(learners, "SuperSurv_grid"),
    identical(unname(fit$event.whichScreen[1L, ]), c(TRUE, FALSE)),
    all(fit$event.whichScreen[2L, ]),
    fit$event.fitLibrary[[1L]]$local_scope_marker == 17L,
    fit$cens.fitLibrary[[1L]]$local_scope_marker == 19L,
    !any(c(fit$event.errorsInCVLibrary, fit$cens.errorsInCVLibrary,
           fit$event.errorsInLibrary, fit$cens.errorsInLibrary)),
    !any(vapply(c(as.character(learners), "surv.local_scope_base",
                   "surv.local_scope_direct", "screen.local_scope"),
                 exists, logical(1L), envir = globalenv(), inherits = FALSE))
  )
  fit
}

check_same <- function(a, b) {
  for (name in c("event.Z", "cens.Z", "event.predict", "cens.predict",
                 "event.coef", "cens.coef", "event.cvRisks", "cens.cvRisks")) {
    stopifnot(isTRUE(all.equal(unname(a[[name]]), unname(b[[name]]), tolerance = 1e-12)))
  }
}

# The same learners through the existing top-level built-in interface.
for (method in c("brier", "logloss")) {
  reference <- do.call(SuperSurv, c(common, list(
    event.library = c("surv.km", "surv.coxph"), cens.library = "surv.km",
    metalearner = method
  )))
  fit <- local_fit(method = method)
  check_same(fit, reference)
  if (method == "brier") {
    brier_fit <- fit
    # A separate verification run installs a multisession future plan before
    # sourcing this file, so this branch also runs on independent workers.
    check_same(local_fit(parallel = TRUE), reference)
  }
}

# The returned object needs saved models and their S3 methods, not local names.
new_X <- X[1:8, ]
new_times <- sort(c(times, as.numeric(quantile(time, c(0.25, 0.45)))))
later <- predict(brier_fit, newdata = new_X, new.times = new_times)
stopifnot(
  max(abs(later$event.predict[1:4, match(times, new_times)] - brier_fit$event.predict)) < 1e-12,
  max(abs(later$cens.predict[1:4, match(times, new_times)] - brier_fit$cens.predict)) < 1e-12
)
input <- tempfile(fileext = ".rds")
output <- tempfile(fileext = ".rds")
saveRDS(list(fit = brier_fit, X = new_X, times = new_times), input)
child <- paste0(
  "args <- commandArgs(TRUE); .libPaths(c(args[3], .libPaths())); ",
  "library(SuperSurv); x <- readRDS(args[1]); ",
  "saveRDS(predict(x$fit, newdata=x$X, new.times=x$times), args[2])"
)
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("-e", shQuote(child), shQuote(input), shQuote(output),
                    shQuote(dirname(find.package("SuperSurv")))))
stopifnot(status == 0L, isTRUE(all.equal(readRDS(output), later, tolerance = 1e-12)))
unlink(c(input, output))

# Caller bindings take precedence over identically named package functions.
masked <- local({
  surv.coxph <- SuperSurv::surv.km
  do.call(SuperSurv, c(common, list(event.library = "surv.coxph", cens.library = "surv.km")))
})
km_reference <- do.call(SuperSurv, c(common, list(
  event.library = "surv.km", cens.library = "surv.km")))
check_same(masked, km_reference)

# A locally generated built-in grid also works as a classed vector directly.
classed <- local({
  km_grid <- create_grid("surv.km", list(dummy = 20260922L))
  do.call(SuperSurv, c(common, list(event.library = km_grid, cens.library = km_grid)))
})
check_same(classed, km_reference)

# Built-ins also resolve when the caller cannot see the attached-package path.
isolated <- list2env(common, parent = baseenv())
namespace_only <- evalq({
  learners <- SuperSurv::create_grid("surv.km", list(dummy = 20260923L))
  SuperSurv::SuperSurv(
    time, event, X, newdata, new.times, event.library = learners,
    cens.library = learners, nFolds = nFolds, cvControl = cvControl,
    control = control
  )
}, isolated)
check_same(namespace_only, km_reference)

error_message <- tryCatch(create_grid("surv.missing_local_base", list(x = 1)),
                          error = conditionMessage)
stopifnot(is.character(error_message),
          grepl("`base_learner` contains an unknown function name", error_message, fixed = TRUE))
message("PASS: local grids, custom bases/screeners, caller precedence, both libraries/losses, future fitting, and fresh-process saved-fit prediction.")
