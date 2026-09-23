library(SuperSurv)

step_curve_at <- getFromNamespace(".step_curve_at", "SuperSurv")
survcomputeCoef <- getFromNamespace(".survcomputeCoef", "SuperSurv")

grid <- c(0, 1, 2)
curve <- c(1, 0.8, 0.5)
stopifnot(
  identical(step_curve_at(grid, curve, c(1, 1.5)), c(0.8, 0.8)),
  identical(
    step_curve_at(grid, curve, c(1, 1.5), left.limit = TRUE),
    c(1, 0.8)
  )
)

preds <- cbind(rep(0, 4), rep(1, 4))
outcome <- c(0.2, 0.4, 0.6, 0.8)
simplex_fit <- getFromNamespace(".simplex_least_squares", "SuperSurv")(
  preds = preds,
  outcome = outcome,
  weights = rep(1, 4),
  tol = 1e-10,
  maxit = 10000L
)
stopifnot(
  simplex_fit$converged,
  all(simplex_fit$weights >= 0),
  abs(sum(simplex_fit$weights) - 1) < 1e-10,
  max(abs(simplex_fit$weights - c(0.5, 0.5))) < 1e-6
)

reference_two_learner_fit <- function(preds, outcome, weights, tol, maxit) {
  weight_sum <- sum(weights)
  gram <- crossprod(preds, weights * preds) / weight_sum
  linear <- drop(crossprod(preds, weights * outcome) / weight_sum)
  largest_eigenvalue <- max(
    eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  )
  current <- c(0.5, 0.5)
  if (!is.finite(largest_eigenvalue) ||
      largest_eigenvalue <= .Machine$double.eps) {
    return(list(
      weights = current,
      converged = TRUE,
      iterations = 0L,
      objective = mean(weights * (outcome - drop(preds %*% current))^2)
    ))
  }

  step <- 1 / (2 * largest_eigenvalue)
  converged <- FALSE
  iterations <- maxit
  for (iteration in seq_len(maxit)) {
    gradient <- 2 * drop(gram %*% current - linear)
    updated <- getFromNamespace(".project_simplex", "SuperSurv")(
      current - step * gradient
    )
    if (max(abs(updated - current)) <= tol) {
      current <- updated
      converged <- TRUE
      iterations <- iteration
      break
    }
    current <- updated
  }
  list(
    weights = current,
    converged = converged,
    iterations = iterations,
    objective = mean(weights * (outcome - drop(preds %*% current))^2)
  )
}

set.seed(93005)
equivalence_predictions <- cbind(
  stats::runif(500, 0.1, 0.9),
  stats::runif(500, 0.1, 0.9)
)
equivalence_outcome <- stats::runif(500, -0.2, 1)
equivalence_weights <- stats::runif(500, 0.5, 1.5)
reference_fit <- reference_two_learner_fit(
  equivalence_predictions, equivalence_outcome, equivalence_weights,
  tol = 1e-8, maxit = 10000L
)
optimized_fit <- getFromNamespace(".simplex_least_squares", "SuperSurv")(
  equivalence_predictions, equivalence_outcome, equivalence_weights,
  tol = 1e-8, maxit = 10000L
)
stopifnot(
  identical(optimized_fit$converged, reference_fit$converged),
  identical(optimized_fit$iterations, reference_fit$iterations),
  max(abs(optimized_fit$weights - reference_fit$weights)) < 1e-12,
  abs(optimized_fit$objective - reference_fit$objective) < 1e-14
)

boundary_predictions <- cbind(rep(0.2, 40), rep(0.8, 40))
boundary_outcome <- rep(0, 40)
reference_boundary_fit <- reference_two_learner_fit(
  boundary_predictions, boundary_outcome, rep(1, 40),
  tol = 1e-8, maxit = 10000L
)
optimized_boundary_fit <- getFromNamespace(
  ".simplex_least_squares", "SuperSurv"
)(
  boundary_predictions, boundary_outcome, rep(1, 40),
  tol = 1e-8, maxit = 10000L
)
stopifnot(
  identical(
    optimized_boundary_fit$iterations,
    reference_boundary_fit$iterations
  ),
  max(abs(
    optimized_boundary_fit$weights - reference_boundary_fit$weights
  )) < 1e-12
)

duplicate_preds <- cbind(outcome, outcome, rev(outcome))
duplicate_fit <- getFromNamespace(".simplex_least_squares", "SuperSurv")(
  preds = duplicate_preds,
  outcome = outcome,
  weights = rep(1, 4),
  tol = 1e-10,
  maxit = 10000L
)
uniform_loss <- mean((outcome - drop(duplicate_preds %*% rep(1 / 3, 3)))^2)
stopifnot(
  duplicate_fit$converged,
  all(duplicate_fit$weights >= 0),
  abs(sum(duplicate_fit$weights) - 1) < 1e-10,
  duplicate_fit$objective <= uniform_loss + 1e-10
)

inclusive <- survcomputeCoef(
  time = 1,
  event = 1,
  t.vals = 1,
  cens.vals = 1,
  preds = matrix(c(0, 1), nrow = 1),
  obsWeights = 1,
  method = "brier",
  event.inclusive = TRUE
)
strict <- survcomputeCoef(
  time = 1,
  event = 1,
  t.vals = 1,
  cens.vals = 1,
  preds = matrix(c(0, 1), nrow = 1),
  obsWeights = 1,
  method = "brier",
  event.inclusive = FALSE
)
stopifnot(inclusive[1] > 0.999, strict[2] > 0.999)

boundary_logloss <- survcomputeCoef(
  time = c(1, 2, 3, 4),
  event = c(1, 0, 1, 0),
  t.vals = rep(2.5, 4),
  cens.vals = c(1e-8, 1, 0.5, 0.5),
  G_t = c(1, 1, 1e-8, 0.5),
  preds = cbind(c(0, 0, 1, 1), c(1, 1, 0, 0)),
  obsWeights = rep(1, 4),
  method = "logloss",
  ipcw.floor = 1e-4,
  ipcw.cap = 100,
  prob.eps = 1e-8
)
boundary_diagnostics <- attr(boundary_logloss, "diagnostics")
stopifnot(
  all(is.finite(boundary_logloss)),
  abs(sum(boundary_logloss) - 1) < 1e-8,
  is.finite(boundary_diagnostics$objective),
  boundary_diagnostics$floor.fraction > 0,
  boundary_diagnostics$cap.fraction > 0,
  identical(boundary_diagnostics$effective.denominator.threshold, 0.01),
  boundary_diagnostics$failure.weights$rows > 0,
  boundary_diagnostics$survivor.weights$rows > 0,
  boundary_diagnostics$combined.weights$cap.count > 0,
  boundary_diagnostics$combined.weights$effective.sample.size > 0,
  boundary_diagnostics$combined.weights$capped.loss.fraction >= 0,
  boundary_diagnostics$combined.weights$capped.loss.fraction <= 1,
  identical(
    boundary_diagnostics$combined.weights$
      stabilized.weight.quantiles[["max"]],
    100
  )
)

control <- SuperSurv:::SuperSurv.control(
  tol = 1e-6,
  ipcw.floor = 1e-3,
  ipcw.cap = 50,
  logloss.eps = 1e-8
)
stopifnot(
  identical(control$tol, 1e-6),
  identical(control$ipcw.floor, 1e-3),
  identical(control$ipcw.cap, 50),
  identical(control$logloss.eps, 1e-8),
  identical(control$time.weighting, "uniform")
)

rmst_inference_error <- tryCatch(
  estimate_marginal_rmst(
    fit = structure(
      list(
        event.coef = 1,
        cens.coef = 1,
        event.libraryNames = data.frame(predAlgorithm = "surv.km", screenAlgorithm = "screen.all"),
        cens.libraryNames = data.frame(predAlgorithm = "surv.km", screenAlgorithm = "screen.all"),
        control = list(),
        varNames = "group"
      ),
      class = "SuperSurv"
    ),
    data = data.frame(group = 0),
    trt_col = "group",
    times = c(1, 2),
    tau = 1,
    inference = TRUE
  ),
  error = function(error) conditionMessage(error)
)
stopifnot(grepl("validated inference procedure", rmst_inference_error, fixed = TRUE))

# Enumerate an independent discrete joint distribution, including T = C ties.
joint <- expand.grid(T = c(1, 2, 4), C = c(1, 3, 5))
joint$probability <- rep(c(0.3, 0.4, 0.3), 3) * rep(c(0.2, 0.3, 0.5), each = 3)
observed <- pmin(joint$T, joint$C)
delta <- as.numeric(joint$T <= joint$C)
S_observed <- vapply(observed, function(t) sum(c(0.3, 0.4, 0.3)[c(1, 2, 4) > t]), numeric(1))
G_left <- vapply(observed, function(t) sum(c(0.2, 0.3, 0.5)[c(1, 3, 5) >= t]), numeric(1))
constant_candidates <- cbind(rep(0, nrow(joint)), rep(1, nrow(joint)))
for (horizon in c(1, 2, 3)) {
  true_S <- sum(c(0.3, 0.4, 0.3)[c(1, 2, 4) > horizon])
  true_G <- sum(c(0.2, 0.3, 0.5)[c(1, 3, 5) > horizon])
  event_fit <- survcomputeCoef(observed, delta, horizon, G_left,
    constant_candidates, joint$probability, ipcw.cap = Inf,
    optimizer.tol = 1e-12)
  censor_fit <- survcomputeCoef(observed, 1 - delta, horizon, S_observed,
    constant_candidates, joint$probability, event.inclusive = TRUE,
    ipcw.cap = Inf, optimizer.tol = 1e-12)
  stopifnot(abs(event_fit[2] - true_S) < 1e-8,
            abs(censor_fit[2] - true_G) < 1e-8)

  # A single candidate must still evaluate the log-loss and its diagnostics.
  probability <- 0.6
  log_fit <- survcomputeCoef(observed, delta, horizon, G_left,
    matrix(probability, nrow(joint), 1), joint$probability,
    method = "logloss", G_t = rep(true_G, nrow(joint)), ipcw.cap = Inf)
  expected <- -(1 - true_S) * log(1 - probability) - true_S * log(probability)
  diagnostic <- attr(log_fit, "diagnostics")
  stopifnot(diagnostic$converged, is.finite(diagnostic$objective),
            abs(diagnostic$objective * nrow(joint) - expected) < 1e-10)
}
