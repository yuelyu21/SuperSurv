library(SuperSurv)

data("metabric", package = "SuperSurv")
dat <- metabric[1:30, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
new_times <- seq(20, 100, by = 20)

expect_error_contains <- function(expr, pattern) {
  msg <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  stopifnot(!is.null(msg), grepl(pattern, msg, fixed = TRUE))
}

km_grid <- create_grid("surv.km", list(dummy = 1))

fit_grid <- SuperSurv(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X[1:5, , drop = FALSE],
  new.times = new_times,
  event.library = km_grid,
  cens.library = "surv.km",
  nFolds = 2,
  verbose = FALSE
)

stopifnot(
  inherits(km_grid, "SuperSurv_grid"),
  is.character(km_grid),
  inherits(fit_grid, "SuperSurv")
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X[1:5, , drop = FALSE],
    new.times = new_times,
    event.library = 1,
    cens.library = "surv.km",
    nFolds = 2,
    verbose = FALSE
  ),
  "`event.library` must be either a character vector of learner names"
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X[1:5, , drop = FALSE],
    new.times = new_times,
    event.library = c("surv.km", ""),
    cens.library = "surv.km",
    nFolds = 2,
    verbose = FALSE
  ),
  "`event.library` contains missing or empty function names"
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X[1:5, , drop = FALSE],
    new.times = new_times,
    event.library = "does.not.exist",
    cens.library = "surv.km",
    nFolds = 2,
    verbose = FALSE
  ),
  "`event.library` contains unknown function name(s): does.not.exist"
)

expect_error_contains(
  SuperSurv(
    time = dat$duration,
    event = dat$event,
    X = X,
    newdata = X[1:5, , drop = FALSE],
    new.times = new_times,
    event.library = "surv.km",
    cens.library = list(list("surv.km")),
    nFolds = 2,
    verbose = FALSE
  ),
  "`cens.library` element 1 must be a character vector"
)
