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
screen.marg <- function(time, event, X, obsWeights, minscreen = 2, min.p = 0.1, ...) {

  requireNamespace("survival", quietly = TRUE)

  pvals <- apply(X, 2, function(col) {
    # Suppress warnings in case a single feature perfectly separates the data
    est <- suppressWarnings(
      survival::coxph(survival::Surv(time, event) ~ col, weights = obsWeights)
    )
    # Extract the Wald test p-value safely
    smry <- summary(est)
    if ("waldtest" %in% names(smry)) {
      return(smry$waldtest['pvalue'])
    } else {
      return(1) # Return high p-value if model failed
    }
  })

  whichVariable <- pvals <= min.p

  # Safety net
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
screen.glmnet <- function(time, event, X, obsWeights, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...) {

  requireNamespace("glmnet", quietly = TRUE)

  if (!is.matrix(X)) {
    X <- stats::model.matrix(~-1 + ., X)
  }

  # Hack to prevent glmnet from crashing on tied times
  time[event == 0] <- time[event == 0] + min(diff(sort(unique(time)))) / 2
  if(any(time == 0)) time[time == 0] <- min(time[time > 0]) / 2

  # Fit CV Lasso
  fit.glmnet <- suppressWarnings(
    glmnet::cv.glmnet(y = survival::Surv(time, event), x = X,
                      weights = obsWeights, family = 'cox',
                      alpha = alpha, nfolds = nfolds, nlambda = nlambda)
  )

  # Extract active variables at optimal lambda
  whichVariable <- (as.numeric(coef(fit.glmnet$glmnet.fit, s = fit.glmnet$lambda.min)) != 0)

  # Safety net
  if (sum(whichVariable) < minscreen) {
    warning("Fewer than minscreen variables passed the glmnet screen, relaxing lambda.")
    sumCoef <- apply(as.matrix(fit.glmnet$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fit.glmnet$glmnet.fit$beta)[,newCut] != 0)
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
screen.var <- function(time, event, X, obsWeights, keep_fraction = 0.5, minscreen = 2, ...) {

  # Calculate variance of each feature
  vars <- apply(X, 2, var, na.rm = TRUE)

  # Determine the cutoff threshold based on the fraction we want to keep
  cutoff <- stats::quantile(vars, probs = 1 - keep_fraction, na.rm = TRUE)

  whichVariable <- vars >= cutoff

  # Safety net
  if(sum(whichVariable) < minscreen) {
    whichVariable <- rep(FALSE, ncol(X))
    whichVariable[order(vars, decreasing = TRUE)[1:minscreen]] <- TRUE
  }

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
