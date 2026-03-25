#' Explain Predictions with Global SHAP (Kernel SHAP)
#'
#' @param model A fitted SuperSurv object OR a single wrapper output.
#' @param X_explain The dataset you want to explain (e.g., `X_test[1:10, ]`).
#' @param X_background The reference dataset for fastshap (e.g., `X_train[1:100, ]`).
#' @param nsim Number of simulations. Defaults to 20.
#' @param only_best Logical. If TRUE and model is SuperSurv, only explains the highest-weighted base learner.
#' @param verbose Logical; if \code{TRUE}, progress messages are shown.
#' @return A data.frame of class \code{c("explain", "data.frame")} containing the calculated SHAP values. The columns correspond to the covariates in \code{X_explain}.
#' @examples
#' if (requireNamespace("fastshap", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   new.times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = new.times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   shap_values <- explain_kernel(
#'     model = fit,
#'     X_explain = X[1:10, , drop = FALSE],
#'     X_background = X[11:40, , drop = FALSE],
#'     nsim = 5
#'   )
#'
#'   dim(shap_values)
#' }
#' @export
explain_kernel <- function(model, X_explain, X_background, nsim = 20,
                           only_best = FALSE,verbose = FALSE) {

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

      if (isTRUE(verbose)) { message(sprintf(" -> SHAP for %s (Weight: %.3f)", model_name, weight))  }

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
    if (isTRUE(verbose)) {  message(  " -> Calculating SHAP for single learner of class: ",  class(model_to_explain)[1]    ) }

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
#' @examples
#' if (requireNamespace("survex", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   times <- seq(20, 120, by = 20)
#'   y <- survival::Surv(dat$duration, dat$event)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   explainer <- explain_survex(
#'     model = fit,
#'     data = X,
#'     y = y,
#'     times = times,
#'     label = "SuperSurv_demo"
#'   )
#'
#'   class(explainer)
#' }
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
    stop("Input must be a fitted 'SuperSurv' object or a valid single learner wrapper output 'surv.*'.")
  }

  # 2. Define the Custom Prediction Wrapper

  custom_predict <- function(mod, newdata, ...) {
    if (inherits(mod, "SuperSurv")) {
      preds <- predict(object = mod, newdata = newdata, new.times = times)
      return(preds$event.predict)
    }

    if (is.list(mod) && !is.null(mod$fit)) {
      return(predict(mod$fit, newdata = newdata, new.times = times))
    }

    stop("Unsupported model type supplied to custom_predict().")
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





#' Plot Global Feature Importance for SuperSurv
#' @param title Plot title.
#' @param shap_values The output from \code{explain_kernel()}.
#' @param top_n Number of features to show (default 10)
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @examples
#' if (requireNamespace("fastshap", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   new.times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = new.times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   shap_values <- explain_kernel(
#'     model = fit,
#'     X_explain = X[1:10, , drop = FALSE],
#'     X_background = X[11:40, , drop = FALSE],
#'     nsim = 5
#'   )
#'
#'   plot_global_importance(shap_values, top_n = 5)
#' }
#' @export
plot_global_importance <- function(shap_values, title = "SuperSurv: Ensemble Feature Importance", top_n = 10) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # 1. Calculate Mean |SHAP|
  # We use drop=FALSE and as.matrix to ensure it works with different data types
  importance <- colMeans(abs(as.matrix(shap_values)))

  # 2. Organize for ggplot
  df_imp <- data.frame(
    Feature = names(importance),
    Importance = as.numeric(importance)
  ) %>%
    dplyr::arrange(desc(Importance)) %>%
    dplyr::slice(1:top_n) # Take top N features

  # 3. Plot
  ggplot2::ggplot(df_imp, ggplot2::aes(x = stats::reorder(Feature, Importance), y = Importance)) +
    ggplot2::geom_col(fill = "steelblue", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = paste("Based on Weighted Ensemble SHAP Values"),
      x = NULL,
      y = "Mean |SHAP| (Impact on Mortality Risk)"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}




#' Beeswarm Summary Plot for SuperSurv SHAP
#'
#' @param shap_values The output from \code{explain_kernel()}.
#' @param data The covariate data used (X_explain)
#' @param top_n Number of features to display
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @examples
#' if (requireNamespace("fastshap", quietly = TRUE) &&
#'     requireNamespace("ggforce", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   new.times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = new.times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   shap_values <- explain_kernel(
#'     model = fit,
#'     X_explain = X[1:20, , drop = FALSE],
#'     X_background = X[21:50, , drop = FALSE],
#'     nsim = 5
#'   )
#'
#'   plot_beeswarm(
#'     shap_values = shap_values,
#'     data = X[1:20, , drop = FALSE],
#'     top_n = 5
#'   )
#' }
#' @export
plot_beeswarm <- function(shap_values, data, top_n = 10) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("ggforce", quietly = TRUE)

  # 1. Convert to standard data frame to avoid tibble issues
  shap_df <- as.data.frame(shap_values)
  data_df <- as.data.frame(data)

  # 2. Calculate overall importance to rank the Y-axis
  # We do this BEFORE pivoting to find the top features
  importance_df <- data.frame(
    Feature = names(shap_df),
    mean_abs = colMeans(abs(shap_df))
  ) %>%
    dplyr::arrange(desc(mean_abs)) %>%
    dplyr::slice(1:min(top_n, ncol(shap_df))) # EXPLICIT CALL

  top_features <- importance_df$Feature

  # 3. Reshape SHAP values to long format (only for top features)
  shap_long <- shap_df[, top_features, drop = FALSE] %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "Feature",
                        values_to = "SHAP")

  # 4. Reshape original data and scale (0 to 1) for color coding
  # We use a small epsilon (1e-9) to avoid division by zero
  data_long <- data_df[, top_features, drop = FALSE] %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                ~ (. - min(., na.rm=T)) / (max(., na.rm=T) - min(., na.rm=T) + 1e-9))) %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "Feature",
                        values_to = "FeatureValue")

  # 5. Merge
  df_plot <- data.frame(
    Feature = shap_long$Feature,
    SHAP = shap_long$SHAP,
    FeatureValue = data_long$FeatureValue
  )

  # Ensure the Y-axis follows the importance ranking
  df_plot$Feature <- factor(df_plot$Feature, levels = rev(top_features))

  # 6. Plot
  ggplot2::ggplot(df_plot, ggplot2::aes(x = SHAP, y = Feature, color = FeatureValue)) +
    ggforce::geom_sina(alpha = 0.7, size = 1.8, method = "counts", maxwidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linetype = "solid", alpha = 0.5) +
    # SHAP standard colors: Blue (Low) to Red (High)
    ggplot2::scale_color_gradient(low = "#1E88E5", high = "#FF0D57",
                                  breaks = c(0, 1), labels = c("Low", "High"),
                                  name = "Feature Value") +
    ggplot2::labs(
      title = "SuperSurv SHAP Summary",
      subtitle = paste("Top", length(top_features), "Features by Mean Absolute SHAP"),
      x = "SHAP Value (Impact on Mortality Risk)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold",size = 16),
      plot.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.title = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 16)
    )

}




#' Waterfall Plot for an Individual Patient
#'
#' @param shap_values The output from \code{explain_kernel()}.
#' @param patient_index The row index of the patient to explain
#' @param top_n Number of features to show (default 10)
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @examples
#' if (requireNamespace("fastshap", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   new.times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = new.times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   shap_values <- explain_kernel(
#'     model = fit,
#'     X_explain = X[1:10, , drop = FALSE],
#'     X_background = X[11:40, , drop = FALSE],
#'     nsim = 5
#'   )
#'
#'   plot_patient_waterfall(
#'     shap_values = shap_values,
#'     patient_index = 1,
#'     top_n = 5
#'   )
#' }
#' @export
plot_patient_waterfall <- function(shap_values, patient_index = 1, top_n = 10) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # 1. Extract the SHAP values for that specific patient
  row_vals <- as.data.frame(shap_values)[patient_index, ]

  df_p <- data.frame(
    Feature = names(row_vals),
    SHAP = as.numeric(row_vals)
  ) %>%
    dplyr::arrange(desc(abs(SHAP))) %>%
    dplyr::slice(1:top_n) %>%
    dplyr::mutate(Direction = ifelse(SHAP > 0, "Higher Risk", "Lower Risk"))

  # 2. Plot
  ggplot2::ggplot(df_p, ggplot2::aes(x = stats::reorder(Feature, SHAP), y = SHAP, fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Higher Risk" = "#FF0D57", "Lower Risk" = "#1E88E5")) +
    ggplot2::labs(
      title = paste("Patient", patient_index, ": Risk Factor Contribution"),
      subtitle = "How individual features pushed this patient away from the average risk",
      x = NULL,
      y = "SHAP Value (Impact on Risk Score)"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}







#' Plot SHAP Dependence for SuperSurv
#'
#' @param shap_values The output from \code{explain_kernel()}.
#' @param data The original covariate data used for the explanation (X_explain)
#' @param feature_name String name of the column to plot
#' @param title Optional custom title.
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @examples
#' if (requireNamespace("fastshap", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   new.times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = new.times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   shap_values <- explain_kernel(
#'     model = fit,
#'     X_explain = X[1:20, , drop = FALSE],
#'     X_background = X[21:50, , drop = FALSE],
#'     nsim = 5
#'   )
#'
#'   plot_dependence(
#'     shap_values = shap_values,
#'     data = X[1:20, , drop = FALSE],
#'     feature_name = colnames(X)[1]
#'   )
#' }
#' @export
plot_dependence <- function(shap_values, data, feature_name, title = NULL) {

  requireNamespace("ggplot2", quietly = TRUE)
  if(is.null(title)) title <- paste("SHAP Dependence:", feature_name)

  df_plot <- data.frame(
    Feature_Value = data[[feature_name]],
    SHAP_Value = shap_values[[feature_name]]
  )

  ggplot2::ggplot(df_plot, ggplot2::aes(x = Feature_Value, y = SHAP_Value)) +
    # 1. Smooth line in background
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "black",
                         se = TRUE, fill = "gray80", linetype = "dashed") +
    # 2. Points on top
    ggplot2::geom_point(alpha = 0.7, color = "firebrick", size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    ggplot2::labs(title = title, x = paste(feature_name, "Value"), y = "SHAP Value") +
    ggplot2::theme_minimal(base_size = 14)
}












#' Survival Probability Heatmap
#'
#' @param object A fitted SuperSurv object
#' @param newdata Test covariates (e.g., `X_te[1:50, ]`)
#' @param times The time grid to visualize
#' @return A \code{ggplot} object visualizing the SHAP values.
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:80, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'   times <- seq(20, 120, by = 20)
#'
#'   fit <- SuperSurv(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     newdata = X,
#'     new.times = times,
#'     event.library = c("surv.coxph", "surv.ridge"),
#'     cens.library = c("surv.coxph"),
#'     control = list(saveFitLibrary = TRUE)
#'   )
#'
#'   plot_survival_heatmap(
#'     object = fit,
#'     newdata = X[1:20, , drop = FALSE],
#'     times = times
#'   )
#' }
#' @export
plot_survival_heatmap <- function(object, newdata, times) {

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # 1. Generate Predictions
  preds <- predict(object, newdata = newdata, new.times = times)$event.predict

  # 2. Prepare Data: Sort patients by their final survival probability
  # This makes the heatmap look "organized" (High risk at top, Low risk at bottom)
  order_idx <- order(preds[, ncol(preds)])
  preds_sorted <- preds[order_idx, ]

  df_heatmap <- as.data.frame(preds_sorted)
  colnames(df_heatmap) <- times
  df_heatmap$PatientID <- factor(1:nrow(df_heatmap))

  # 3. Pivot to Long Format
  df_long <- df_heatmap %>%
    tidyr::pivot_longer(cols = -PatientID, names_to = "Time", values_to = "Survival") %>%
    dplyr::mutate(Time = as.numeric(Time))

  ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = PatientID, fill = Survival)) +
    ggplot2::geom_tile() +
    # Using the 'Plasma' palette to match the purple/red/yellow theme of survex
    ggplot2::scale_fill_viridis_c(option = "plasma", direction = 1, name = "S(t)") +
    ggplot2::labs(
      title = "SuperSurv: Patient Risk Trajectories",
      subtitle = "Sorted by predicted survival at final follow-up",
      x = "Time (Months)",
      y = "Individual Patients"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      # Match the dark purple text from your importance plots
      text = ggplot2::element_text(color = "#371ea3"),
      plot.title = ggplot2::element_text(face = "bold", size = 18)
    )
}
