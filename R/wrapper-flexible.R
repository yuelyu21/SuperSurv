.check_wrapper_dots <- function(dots, reserved, learner) {
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                       anyDuplicated(names(dots)))) {
    stop("Additional arguments for `", learner,
         "` must have unique, non-empty names.", call. = FALSE)
  }
  conflict <- intersect(names(dots), reserved)
  if (length(conflict)) {
    stop("`", learner, "` manages or does not support these arguments: ",
         paste(conflict, collapse = ", "), ".", call. = FALSE)
  }
  invisible(dots)
}

.wrapper_feature_spec <- function(X, learner) {
  if (any(names(X) %in% c(".supersurv_time", ".supersurv_event"))) {
    stop("Feature names `.supersurv_time` and `.supersurv_event` are reserved by `",
         learner, "`.", call. = FALSE)
  }
  list(features = names(X),
       factor_levels = lapply(X[vapply(X, is.factor, logical(1L))], levels))
}

.wrapper_prediction_data <- function(newdata, spec, learner) {
  if (!is.data.frame(newdata) || !setequal(names(newdata), spec$features)) {
    stop("`newdata` must be a data frame containing exactly the features used by `",
         learner, "`: ", paste(spec$features, collapse = ", "), ".",
         call. = FALSE)
  }
  .validate_data_frame_columns(newdata, "newdata")
  newdata <- newdata[, spec$features, drop = FALSE]
  for (feature in names(spec$factor_levels)) {
    unknown <- setdiff(as.character(newdata[[feature]]), spec$factor_levels[[feature]])
    if (length(unknown)) {
      stop("`newdata` contains unseen level(s) for `", feature, "` in `",
           learner, "`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    newdata[[feature]] <- factor(newdata[[feature]],
                                levels = spec$factor_levels[[feature]],
                                ordered = is.ordered(newdata[[feature]]))
  }
  newdata
}

#' Flexible Parametric Distribution Survival Learner
#'
#' Fits a weighted parametric survival model using [flexsurv::flexsurvreg()].
#' Native survival probabilities are used without risk-score calibration.
#'
#' @inheritParams surv.mboost
#' @param dist Distribution: `"gengamma"` (generalized gamma, default),
#'   `"gompertz"`, `"gamma"`, `"weibull"`, `"exp"`, `"lnorm"`, or `"llogis"`.
#' @param ... Additional named arguments to [flexsurv::flexsurvreg()], such as
#'   `inits`, `anc`, or optimizer `control`. Outcome, data, weight, truncation,
#'   and relative-survival arguments are managed or excluded by this adapter.
#' @details Only single-event, right-censored outcomes are supported. Covariates
#'   enter the distribution's location parameter by default; `anc` can specify
#'   covariates on ancillary parameters. Generalized gamma may require suitable
#'   initial values, especially in small training folds. Failed optimization is
#'   reported as an error rather than silently accepted.
#' @return A list with numeric survival matrix `pred` and fitted object `fit`.
#' @examples
#' if (requireNamespace("flexsurv", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   X <- dat[, "x1", drop = FALSE]
#'   fit <- surv.flexsurvreg(dat$duration, dat$event, X,
#'                          X[1:3, , drop = FALSE], c(50, 100), dist = "weibull")
#'   predict(fit$fit, X[1:3, , drop = FALSE], new.times = c(25, 50, 100))
#' }
#' @export
surv.flexsurvreg <- function(time, event, X, newdata = NULL, new.times,
                             obsWeights = NULL, id = NULL,
                             dist = "gengamma", ...) {
  learner <- "surv.flexsurvreg"
  input <- .prepare_wrapper_inputs(time, event, X, newdata, new.times,
                                   obsWeights, learner)
  choices <- c("gengamma", "gompertz", "gamma", "weibull", "exp", "lnorm", "llogis")
  if (!is.character(dist) || length(dist) != 1L || is.na(dist) ||
      !dist %in% choices) {
    stop("`dist` must be one of: ", paste(choices, collapse = ", "), ".",
         call. = FALSE)
  }
  spec <- .wrapper_feature_spec(input$X, learner)
  dots <- list(...)
  .check_wrapper_dots(dots, c("formula", "data", "weights", "bhazard", "rtrunc",
                             "subset", "na.action"), learner)
  if (!requireNamespace("flexsurv", quietly = TRUE)) {
    stop("Learner `surv.flexsurvreg` requires the optional package 'flexsurv'.",
         call. = FALSE)
  }
  training <- data.frame(.supersurv_time = input$time,
                         .supersurv_event = input$event, input$X, check.names = FALSE)
  formula <- stats::as.formula("survival::Surv(.supersurv_time, .supersurv_event) ~ .")
  fit <- do.call(flexsurv::flexsurvreg,
                 c(list(formula = formula, data = training,
                        weights = input$obsWeights, dist = dist), dots))
  if (!is.null(fit$opt$convergence) && fit$opt$convergence != 0L) {
    stop("`surv.flexsurvreg` did not converge; inspect `dist`, initial values, ",
         "and optimizer controls.", call. = FALSE)
  }
  fit_object <- list(object = fit, feature_spec = spec, dist = dist)
  class(fit_object) <- learner
  list(pred = predict.surv.flexsurvreg(fit_object, input$newdata, input$new.times),
       fit = fit_object)
}

#' @noRd
#' @export
predict.surv.flexsurvreg <- function(object, newdata, new.times, ...) {
  newdata <- .wrapper_prediction_data(newdata, object$feature_spec, "surv.flexsurvreg")
  new.times <- .validate_time_grid(new.times, "new.times")
  .flexsurv_prediction_matrix(object$object, newdata, new.times, "surv.flexsurvreg")
}

#' Penalized Smooth Hazard Survival Learner
#'
#' Fits an overall hazard model using [survPen::survPen()] and returns its native
#' survival probabilities. Only single-event, right-censored outcomes and uniform
#' observation weights are supported.
#'
#' @inheritParams surv.mboost
#' @param formula Optional one-sided hazard formula using feature names and
#'   `.supersurv_time`. For example,
#'   `~ smf(.supersurv_time, df = 4) + smf(x1, df = 4) + x2`, or
#'   `~ tensor(.supersurv_time, x1, df = c(4, 4)) + x2` for a time-varying effect.
#'   The `survPen` constructors `smf`, `tensor`, `tint`, and `rd` are available
#'   without attaching the backend package. Formula variables must come from
#'   the current training features or `.supersurv_time`.
#' @param baseline.df Baseline smooth degrees of freedom, at least 3. Used only
#'   when `formula` is NULL; the default adds linear terms for all features to
#'   `smf(.supersurv_time, df = baseline.df)`.
#' @param n.legendre Positive integer quadrature order for both fitting and
#'   survival prediction. Increase this to check integration accuracy.
#' @param ... Additional named fitting controls passed to [survPen::survPen()],
#'   such as `lambda`, `method`, or `max.it.beta`.
#' @details Nonuniform observation weights are rejected because the backend does
#'   not expose case-weighted fitting. The adapter does not implement net survival,
#'   relative mortality, or left truncation. Survival curves that increase beyond
#'   numerical tolerance cause an error; increase quadrature accuracy rather than
#'   silently projecting a materially invalid curve. Custom formulas must remain
#'   valid after feature screening; use `screen.all` when explicitly naming features.
#' @return A list with numeric survival matrix `pred` and fitted object `fit`.
#' @examples
#' if (requireNamespace("survPen", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   X <- dat[, "x1", drop = FALSE]
#'   fit <- surv.survPen(dat$duration, dat$event, X,
#'                       X[1:3, , drop = FALSE], c(50, 100), baseline.df = 3)
#'   dim(fit$pred)
#' }
#' @export
surv.survPen <- function(time, event, X, newdata = NULL, new.times,
                         obsWeights = NULL, id = NULL, formula = NULL,
                         baseline.df = 4L, n.legendre = 50L, ...) {
  learner <- "surv.survPen"
  input <- .prepare_wrapper_inputs(time, event, X, newdata, new.times,
                                   obsWeights, learner)
  if (length(unique(input$obsWeights)) > 1L) {
    stop("`surv.survPen` does not support nonuniform `obsWeights`; ",
         "the 'survPen' backend does not expose case-weighted fitting.", call. = FALSE)
  }
  for (argument in c("baseline.df", "n.legendre")) {
    value <- get(argument)
    minimum <- if (argument == "baseline.df") 3L else 1L
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < minimum || value > .Machine$integer.max || value != floor(value)) {
      stop("`", argument, "` must be one integer at least ", minimum, ".", call. = FALSE)
    }
  }
  spec <- .wrapper_feature_spec(input$X, learner)
  dots <- list(...)
  .check_wrapper_dots(dots, c("data", "t1", "t0", "event", "expected", "type",
                             "weights", "cluster"), learner)
  if (!is.null(formula) && (!inherits(formula, "formula") || length(formula) != 2L)) {
    stop("`formula` must be NULL or a one-sided formula for the log hazard.", call. = FALSE)
  }
  if (!requireNamespace("survPen", quietly = TRUE)) {
    stop("Learner `surv.survPen` requires the optional package 'survPen'.", call. = FALSE)
  }
  if (is.null(formula)) {
    terms <- c(paste0("smf(.supersurv_time, df = ", as.integer(baseline.df), ")"),
               vapply(names(input$X), function(x) paste(deparse(as.name(x), backtick = TRUE), collapse = ""), ""))
    formula <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  }
  unknown <- setdiff(all.vars(formula), c(names(input$X), ".supersurv_time"))
  if (length(unknown)) {
    stop("`formula` contains variables not in the training features or ",
         "`.supersurv_time`: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  # survPen parses its smooth terms in its own namespace, including tensor/tint.
  training <- data.frame(.supersurv_time = input$time,
                         .supersurv_event = input$event, input$X, check.names = FALSE)
  fit <- do.call(survPen::survPen,
                 c(list(formula = formula, data = training,
                        t1 = quote(.supersurv_time), event = quote(.supersurv_event),
                        n.legendre = as.integer(n.legendre)), dots))
  if (isFALSE(fit$converged)) {
    stop("`surv.survPen` did not converge; simplify the formula or inspect fitting controls.",
         call. = FALSE)
  }
  fit_object <- list(object = fit, feature_spec = spec, n.legendre = as.integer(n.legendre))
  class(fit_object) <- learner
  list(pred = predict.surv.survPen(fit_object, input$newdata, input$new.times),
       fit = fit_object)
}

#' @noRd
#' @export
predict.surv.survPen <- function(object, newdata, new.times,
                                n.legendre = object$n.legendre, ...) {
  newdata <- .wrapper_prediction_data(newdata, object$feature_spec, "surv.survPen")
  new.times <- .validate_time_grid(new.times, "new.times")
  if (!is.numeric(n.legendre) || length(n.legendre) != 1L ||
      !is.finite(n.legendre) || n.legendre < 1 ||
      n.legendre > .Machine$integer.max || n.legendre != floor(n.legendre)) {
    stop("`n.legendre` must be one positive integer.", call. = FALSE)
  }
  prediction <- matrix(1, nrow(newdata), length(new.times))
  for (j in which(new.times > 0)) {
    data <- newdata
    data$.supersurv_time <- new.times[j]
    prediction[, j] <- stats::predict(object$object, newdata = data,
                                      n.legendre = as.integer(n.legendre),
                                      do.surv = TRUE)$surv
  }
  if (ncol(prediction) > 1L && any(prediction[, -1L, drop = FALSE] -
      prediction[, -ncol(prediction), drop = FALSE] > 1e-8, na.rm = TRUE)) {
    stop("`surv.survPen` returned increasing survival probabilities; ",
         "increase `n.legendre` and check the fitted model.", call. = FALSE)
  }
  .finalize_wrapper_survival(prediction, nrow(newdata), new.times, "surv.survPen")
}
