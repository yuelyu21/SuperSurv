## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----custom-learner-----------------------------------------------------------
my.surv.km <- function(time, event, X, newdata, new.times,
                       obsWeights = NULL, id = NULL, ...) {
  if (is.null(obsWeights)) {
    obsWeights <- rep(1, length(time))
  }

  fit_km <- survival::survfit(
    survival::Surv(time, event) ~ 1,
    weights = obsWeights
  )

  step_surv <- stats::stepfun(fit_km$time, c(1, fit_km$surv), right = FALSE)
  surv_probs <- step_surv(new.times)

  pred <- matrix(
    surv_probs,
    nrow = nrow(newdata),
    ncol = length(new.times),
    byrow = TRUE
  )

  fit <- list(object = fit_km)
  class(fit) <- "my.surv.km"

  list(pred = pred, fit = fit)
}

predict.my.surv.km <- function(object, newdata, new.times, ...) {
  fit_km <- object$object
  step_surv <- stats::stepfun(fit_km$time, c(1, fit_km$surv), right = FALSE)
  surv_probs <- step_surv(new.times)

  matrix(
    surv_probs,
    nrow = nrow(newdata),
    ncol = length(new.times),
    byrow = TRUE
  )
}

## ----custom-learner-check-----------------------------------------------------
data("metabric", package = "SuperSurv")
dat <- metabric[1:40, ]
x_cols <- grep("^x", names(dat))[1:3]
X <- dat[, x_cols, drop = FALSE]
times <- seq(50, 150, by = 50)

wrapper_out <- my.surv.km(
  time = dat$duration,
  event = dat$event,
  X = X,
  newdata = X[1:5, , drop = FALSE],
  new.times = times
)

dim(wrapper_out[["pred"]])
dim(predict(wrapper_out[["fit"]], newdata = X[1:5, , drop = FALSE], new.times = times))

## ----custom-fit, eval=FALSE---------------------------------------------------
# fit <- SuperSurv(
#   time = train$duration,
#   event = train$event,
#   X = X_train,
#   newdata = X_test,
#   new.times = seq(50, 300, by = 50),
#   event.library = c("surv.coxph", "my.surv.km"),
#   cens.library = c("surv.coxph"),
#   control = list(saveFitLibrary = TRUE)
# )
# 
# event_weights(fit)

## ----custom-screener, eval=FALSE----------------------------------------------
# fit <- SuperSurv(
#   time = train$duration,
#   event = train$event,
#   X = X_train,
#   newdata = X_test,
#   new.times = seq(50, 300, by = 50),
#   event.library = list(c("surv.coxph", "screen.all", "my.screen.first3")),
#   cens.library = c("surv.coxph")
# )
# 
# selected_variables(fit, learner = 2)

