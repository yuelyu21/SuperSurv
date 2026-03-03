#' Keep All Variables Screener
#' @param X Training covariate data.frame.
#' @param ... Additional ignored arguments.
#'
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @export
screen.all <- function(X, ...) {
  rep(TRUE, ncol(X))
}




#' Marginal Cox Regression Screening
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param obsWeights Observation weights.
#' @param minscreen Minimum number of covariates to return. Defaults to 2.
#' @param min.p Threshold p-value. Defaults to 0.1.
#' @param ... Additional ignored arguments.
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @export
screen.marg <- function(time, event, X, obsWeights = NULL, minscreen = 2, min.p = 0.1, ...) {

  requireNamespace("survival", quietly = TRUE)

  # vapply safely iterates over data.frame columns without converting them to text
  pvals <- vapply(X, function(col) {

    # Catch hard mathematical errors, not just warnings
    fit_try <- try({
      if (is.null(obsWeights)) {
        survival::coxph(survival::Surv(time, event) ~ col)
      } else {
        survival::coxph(survival::Surv(time, event) ~ col, weights = obsWeights)
      }
    }, silent = TRUE)

    # If the model crashed (e.g., constant variable), assign a p-value of 1 (drop it)
    if (inherits(fit_try, "try-error")) return(1.0)

    # Extract the Wald test p-value safely
    smry <- summary(fit_try)
    if ("waldtest" %in% names(smry)) {
      return(as.numeric(smry$waldtest['pvalue']))
    } else {
      return(1.0)
    }
  }, numeric(1)) # Ensure the output is strictly numeric

  whichVariable <- (pvals <= min.p)

  # Safety net: If everything was dropped, at least keep the top 'minscreen' variables
  if(sum(whichVariable, na.rm = TRUE) < minscreen) {
    whichVariable <- rep(FALSE, ncol(X))
    whichVariable[order(pvals)[1:minscreen]] <- TRUE
  }

  return(whichVariable)
}




#' GLMNET (Lasso) Screening
#' @param time Observed follow-up time.
#' @param event Observed event indicator.
#' @param X Training covariate data.frame.
#' @param obsWeights Observation weights.
#' @param alpha Penalty exponent (1 = lasso).
#' @param minscreen Minimum number of covariates to return. Defaults to 2.
#' @param nfolds Number of CV folds.
#' @param nlambda Number of penalty parameters.
#' @param ... Additional ignored arguments.
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @export
screen.glmnet <- function(time, event, X, obsWeights = NULL, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...) {

  requireNamespace("glmnet", quietly = TRUE)

  # Safe matrix conversion without adding dummy columns
  if (!is.matrix(X)) {
    X <- as.matrix(sapply(X, as.numeric))
  }

  # Safe weights handling
  if (is.null(obsWeights)) {
    obsWeights <- rep(1, nrow(X))
  }

  # Hack to prevent glmnet from crashing on tied times
  if (any(event == 0)) {
    time_diffs <- diff(sort(unique(time)))
    if (length(time_diffs) > 0) {
      time[event == 0] <- time[event == 0] + min(time_diffs) / 2
    }
  }
  if(any(time == 0)) time[time == 0] <- min(time[time > 0]) / 2

  # Fit CV Lasso
  fit.glmnet <- suppressWarnings(
    glmnet::cv.glmnet(y = survival::Surv(time, event), x = X,
                      weights = obsWeights, family = 'cox',
                      alpha = alpha, nfolds = nfolds, nlambda = nlambda)
  )

  # Extract active variables correctly
  coefs <- as.numeric(coef(fit.glmnet, s = "lambda.min"))
  whichVariable <- (coefs != 0)

  # Safety net
  if (sum(whichVariable) < minscreen) {
    sumCoef <- apply(as.matrix(fit.glmnet$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    if (length(newCut) == 0 || is.na(newCut)) {
      whichVariable <- rep(TRUE, ncol(X)) # Ultimate fallback
    } else {
      whichVariable <- (as.matrix(fit.glmnet$glmnet.fit$beta)[,newCut] != 0)
    }
  }

  return(whichVariable)
}




#' Random Survival Forest Screening Algorithm
#'
#' This screening algorithm uses the \code{randomForestSRC} package to select covariates
#' based on their Variable Importance (VIMP). It grows a fast forest and retains features
#' with a VIMP greater than zero.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param X Training covariate data.frame or matrix.
#' @param obsWeights Numeric vector of observation weights.
#' @param minscreen Integer. Minimum number of covariates to return. Defaults to 2.
#' @param ntree Integer. Number of trees to grow. Defaults to 100 for fast screening.
#' @param ... Additional arguments passed to \code{\link[randomForestSRC]{rfsrc}}.
#'
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @export
screen.rfsrc <- function(time, event, X, obsWeights, minscreen = 2, ntree = 100, ...) {
  requireNamespace("randomForestSRC", quietly = TRUE)

  # Fit a fast forest for screening
  data_rf <- data.frame(time = time, event = event, X)
  fit_rf <- randomForestSRC::rfsrc(survival::Surv(time, event) ~ .,
                                   data = data_rf,
                                   ntree = ntree,
                                   case.wt = obsWeights,
                                   importance = TRUE,
                                   nsplit = 2) # Fast splitting

  # Extract Variable Importance (VIMP)
  vimp <- fit_rf$importance

  # Keep variables with positive importance
  whichVariable <- vimp > 0

  # Safety net
  if(sum(whichVariable) < minscreen) {
    whichVariable <- rep(FALSE, ncol(X))
    whichVariable[order(vimp, decreasing = TRUE)[1:minscreen]] <- TRUE
  }

  return(whichVariable)
}



#' High Variance Screening Algorithm (Unsupervised)
#'
#' An unsupervised screening algorithm that filters out low-variance features.
#' This is particularly useful for high-dimensional genomic or transcriptomic data
#' where many features remain relatively constant across all observations.
#'
#' @param time Numeric vector of observed follow-up times (Ignored internally).
#' @param event Numeric vector of event indicators (Ignored internally).
#' @param X Training covariate data.frame or matrix.
#' @param obsWeights Numeric vector of observation weights (Ignored internally).
#' @param keep_fraction Numeric value between 0 and 1. The fraction of highest-variance
#' features to retain. Defaults to 0.5 (keeps the top 50%).
#' @param minscreen Integer. Minimum number of covariates to return. Defaults to 2.
#' @param ... Additional ignored arguments.
#'
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @export
screen.var <- function(time, event, X, obsWeights = NULL, keep_fraction = 0.5, minscreen = 2, ...) {

  # Safely calculate variance of each feature, forcing numeric conversion
  vars <- vapply(X, function(col) var(as.numeric(col), na.rm = TRUE), numeric(1))

  # Determine the cutoff threshold based on the fraction we want to keep
  cutoff <- stats::quantile(vars, probs = 1 - keep_fraction, na.rm = TRUE)
  whichVariable <- (vars >= cutoff)

  # Safety net
  if(sum(whichVariable, na.rm = TRUE) < minscreen) {
    whichVariable <- rep(FALSE, ncol(X))
    whichVariable[order(vars, decreasing = TRUE)[1:minscreen]] <- TRUE
  }

  # Replace any NAs created by zero-variance or all-NA columns with FALSE
  whichVariable[is.na(whichVariable)] <- FALSE

  return(whichVariable)
}





#' Elastic Net Screening Algorithm
#'
#' This screening algorithm uses \code{\link[glmnet]{cv.glmnet}} to select covariates.
#' Unlike LASSO (\code{alpha = 1}), which drops correlated features, Elastic Net
#' (\code{alpha = 0.5} by default) shrinks correlated groups of features together,
#' making it ideal for selecting entire biological pathways.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (1 = event, 0 = censored).
#' @param X Training covariate data.frame or matrix.
#' @param obsWeights Numeric vector of observation weights.
#' @param alpha Numeric penalty exponent for \code{glmnet}. Defaults to 0.5 (Elastic Net).
#' @param minscreen Integer. Minimum number of covariates to return. Defaults to 2.
#' @param nfolds Integer. Number of folds for cross-validation. Defaults to 10.
#' @param nlambda Integer. Number of penalty parameters to search over. Defaults to 100.
#' @param ... Additional arguments passed to \code{screen.glmnet}.
#'
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @export
screen.elasticnet <- function(time, event, X, obsWeights, alpha = 0.5, minscreen = 2, nfolds = 10, nlambda = 100, ...) {

  # We simply pass the arguments directly to our robust screen.glmnet wrapper,
  # forcing alpha = 0.5 instead of alpha = 1.

  whichVariable <- screen.glmnet(
    time = time,
    event = event,
    X = X,
    obsWeights = obsWeights,
    alpha = alpha,
    minscreen = minscreen,
    nfolds = nfolds,
    nlambda = nlambda,
    ...
  )

  return(whichVariable)
}
