#' Explain Predictions with Global SHAP (Kernel SHAP)
#'
#' @param model A fitted SuperSurv object OR a single wrapper output.
#' @param X_explain The dataset you want to explain (e.g., `X_test[1:10, ]`).
#' @param X_background The reference dataset for fastshap (e.g., `X_train[1:100, ]`).
#' @param nsim Number of simulations. Defaults to 20.
#' @param only_best Logical. If TRUE and model is SuperSurv, only explains the highest-weighted base learner.
#' @return A data.frame of class \code{c("explain", "data.frame")} containing the calculated SHAP values. The columns correspond to the covariates in \code{X_explain}.
#' @export
explain_shap <- function(model, X_explain, X_background, nsim = 20, only_best = FALSE) {

  requireNamespace("fastshap", quietly = TRUE)

  # ------------------------------------------------------------------
  # PATH A: The SuperSurv Ensemble
  # ------------------------------------------------------------------
  if (inherits(model, "SuperSurv")) {

    if (only_best) {
      active_indices <- which.max(model$event.coef)
      active_weights <- 1
    } else {
      active_indices <- which(model$event.coef > 0.001)
      active_weights <- model$event.coef[active_indices]
      active_weights <- active_weights / sum(active_weights) # Re-normalize
    }

    shap_list <- list()
    for (i in seq_along(active_indices)) {
      idx <- active_indices[i]
      model_fit <- model$event.fitLibrary[[idx]]
      weight <- active_weights[i]
      model_name <- names(model$event.fitLibrary)[idx]

      message(sprintf(" -> SHAP for %s (Weight: %.3f)", model_name, weight))

      s <- fastshap::explain(
        object = model_fit,
        X = X_explain,
        nsim = nsim,
        adjust = TRUE,
        baseline = mean(get_risk_universal(model_fit, X_background)),
        pred_wrapper = function(obj, newdata) { get_risk_universal(obj, newdata) }
      )
      shap_list[[i]] <- as.data.frame(s) * weight
    }

    final_shap <- Reduce("+", shap_list)
    class(final_shap) <- c("explain", "data.frame")
    return(final_shap)

    # ------------------------------------------------------------------
    # PATH B: A Single Base Learner
    # ------------------------------------------------------------------
  } else if (is.list(model) && !is.null(model$fit)) {

    model_to_explain <- model$fit
    message(" -> Calculating SHAP for single learner of class: ", class(model_to_explain)[1])

    s <- fastshap::explain(
      object = model_to_explain,
      X = X_explain,
      nsim = nsim,
      adjust = TRUE,
      baseline = mean(get_risk_universal(model_to_explain, X_background)),
      pred_wrapper = function(obj, newdata) { get_risk_universal(obj, newdata) }
    )

    shap_out <- as.data.frame(s)
    class(shap_out) <- c("explain", "data.frame")
    return(shap_out)

  } else {
    stop("Input must be a fitted 'SuperSurv' object or a valid single learner wrapper output.")
  }
}




#' Create a Time-Dependent Survex Explainer
#'
#' Bridges a fitted SuperSurv ensemble or a single base learner to the
#' \code{survex} package for Time-Dependent SHAP and Model Parts.
#'
#' @param model A fitted SuperSurv object OR a single wrapper output.
#' @param data Covariate data for explanation (data.frame).
#' @param y The survival object (\code{Surv(time, event)}).
#' @param times The time grid for evaluation.
#' @param label Optional character string to name the explainer.
#' @return An explainer object of class \code{survex_explainer} created by
#'   \code{\link[survex]{explain_survival}}, which can be passed to DALEX and survex functions for further model diagnostics and plotting.
#' @export
explain_survex <- function(model, data, y, times, label = NULL) {

  if (!requireNamespace("survex", quietly = TRUE)) stop("Please install the 'survex' package.")

  # 1. Determine Model & Label
  if (inherits(model, "SuperSurv")) {
    model_obj <- model
    if (is.null(label)) label <- "SuperSurv_Ensemble"
  } else if (is.list(model) && !is.null(model$fit)) {
    model_obj <- model
    if (is.null(label)) label <- class(model$fit)[1]
  } else {
    stop("Input must be a fitted 'SuperSurv' object or a valid single learner wrapper output.")
  }

  # 2. Define the Custom Prediction Wrapper
  custom_predict <- function(mod, newdata, ...) {
    if (inherits(mod, "SuperSurv")) {
      preds <- predict(object = mod, newdata = newdata, new.times = times)
      return(preds$event.SL.predict)
    } else {
      p <- predict(mod$fit, newX = newdata, new.times = times)
      return(p)
    }
  }

  # 3. Create the Explainer
  explainer <- survex::explain_survival(
    model = model_obj,
    data = as.data.frame(data),
    y = y,
    predict_survival_function = custom_predict,
    times = times,
    label = label,
    verbose = FALSE
  )

  return(explainer)
}
