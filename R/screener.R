#' @keywords internal
#' @noRd
.make_design_matrix <- function(X) {
  X_df <- as.data.frame(X)

  # If user passed a matrix, preserve colnames
  if (is.matrix(X)) X_df <- as.data.frame(X)

  # model.matrix handles factors/characters safely via one-hot encoding
  mm <- stats::model.matrix(~ . - 1, data = X_df)

  # Ensure it's a plain numeric matrix
  mm <- as.matrix(mm)
  storage.mode(mm) <- "double"
  mm
}



#' Keep All Variables Screener
#' @param X Training covariate data.frame.
#' @param ... Additional ignored arguments.
#'
#' @return A logical vector of the same length as the number of columns in \code{X},
#'   indicating which variables passed the screening algorithm (\code{TRUE} to keep,
#'   \code{FALSE} to drop).
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:20, ]
#' x_cols <- grep("^x", names(dat))[1:5]
#' X <- dat[, x_cols, drop = FALSE]
#'
#' screen.all(X)
#' @export
screen.all <- function(X, ...) {
  whichVariable <- rep(TRUE, ncol(X))
  names(whichVariable) <- colnames(X)
  return(whichVariable)
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
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:40, ]
#' x_cols <- grep("^x", names(dat))[1:5]
#' X <- dat[, x_cols, drop = FALSE]
#'
#' screen.marg(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   minscreen = 2,
#'   min.p = 0.2
#' )
#' @export
screen.marg <- function(time, event, X, obsWeights = NULL,
                        minscreen = 2, min.p = 0.1, ...) {

  requireNamespace("survival", quietly = TRUE)

  X_mat <- .make_design_matrix(X)

  pvals <- vapply(seq_len(ncol(X_mat)), function(j) {
    colj <- X_mat[, j]

    fit_try <- try({
      if (is.null(obsWeights)) {
        survival::coxph(survival::Surv(time, event) ~ colj)
      } else {
        survival::coxph(survival::Surv(time, event) ~ colj, weights = obsWeights)
      }
    }, silent = TRUE)

    if (inherits(fit_try, "try-error")) return(1.0)

    smry <- summary(fit_try)
    # Single coefficient p-value
    p <- suppressWarnings(as.numeric(smry$coefficients[,"Pr(>|z|)"][1]))
    if (is.na(p)) p <- 1.0
    p
  }, numeric(1))

  whichVariable <- (pvals <= min.p)

  # Safety net: keep at least minscreen
  if (sum(whichVariable, na.rm = TRUE) < minscreen) {
    whichVariable <- rep(FALSE, ncol(X_mat))
    whichVariable[order(pvals)[seq_len(minscreen)]] <- TRUE
  }
  names(whichVariable) <- colnames(X_mat)

  # IMPORTANT: return mask in original X_mat space (design matrix columns)
  # This is OK if downstream uses the same design matrix builder.
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
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:40, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'
#'   screen.glmnet(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     alpha = 1,
#'     minscreen = 2,
#'     nfolds = 3,
#'     nlambda = 20
#'   )
#' }
#' @export
screen.glmnet <- function(time, event, X, obsWeights = NULL,
                          alpha = 1, minscreen = 2,
                          nfolds = 10, nlambda = 100, ...) {

  requireNamespace("glmnet", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)

  # ---- Fix 1: Cox family requires strictly positive times ----
  time <- pmax(time, 1e-5)

  # Build a safe numeric design matrix (same philosophy as your wrappers)
  X_df <- as.data.frame(X)
  X_mat <- stats::model.matrix(~ . - 1, data = X_df)
  X_mat <- as.matrix(X_mat)
  storage.mode(X_mat) <- "double"

  if (is.null(obsWeights)) obsWeights <- rep(1, nrow(X_mat))

  fit.glmnet <- suppressWarnings(
    glmnet::cv.glmnet(
      x = X_mat,
      y = survival::Surv(time, event),
      weights = obsWeights,
      family = "cox",
      alpha = alpha,
      nfolds = nfolds,
      nlambda = nlambda
    )
  )

  # ---- Fix 2: use exported coef() ----
  beta_min <- as.matrix(stats::coef(fit.glmnet, s = "lambda.min"))
  whichVariable <- as.vector(beta_min != 0)

  # Safety net: keep at least minscreen
  if (sum(whichVariable) < minscreen) {
    beta_path <- as.matrix(fit.glmnet$glmnet.fit$beta)  # p x nlambda
    nnz <- apply(beta_path, 2, function(b) sum(b != 0))

    idx <- which(nnz >= minscreen)
    if (length(idx) > 0) {
      use <- idx[1]
      whichVariable <- as.vector(beta_path[, use] != 0)
    } else {
      use <- which.max(nnz)
      whichVariable <- as.vector(beta_path[, use] != 0)
      if (sum(whichVariable) == 0) whichVariable <- rep(TRUE, ncol(X_mat))
    }
  }

  names(whichVariable) <- colnames(X_mat)

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
#' @examples
#' if (requireNamespace("randomForestSRC", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:40, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'
#'   screen.rfsrc(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     minscreen = 2,
#'     ntree = 10
#'   )
#' }
#' @export
screen.rfsrc <- function(time, event, X, obsWeights = NULL,
                         minscreen = 2, ntree = 100, ...) {

  requireNamespace("randomForestSRC", quietly = TRUE)
  requireNamespace("survival", quietly = TRUE)
  Surv <- survival::Surv

  X_df <- as.data.frame(X)

  if (is.null(obsWeights)) obsWeights <- rep(1, nrow(X_df))

  data_rf <- data.frame(time = time, event = event, X_df)

  fit_rf <- randomForestSRC::rfsrc(
    Surv(time, event) ~ .,
    data = data_rf,
    ntree = ntree,
    case.wt = obsWeights,
    importance = TRUE,
    nsplit = 2,
    ...
  )

  vimp <- fit_rf$importance
  if (is.null(vimp)) {
    # fallback: keep all if VIMP is unavailable
    return(rep(TRUE, ncol(X_df)))
  }

  # vimp is named by original X columns (not model.matrix columns)
  vimp <- as.numeric(vimp)
  whichVariable <- vimp > 0
  whichVariable[is.na(whichVariable)] <- FALSE

  if (sum(whichVariable) < minscreen) {
    whichVariable <- rep(FALSE, length(vimp))
    whichVariable[order(vimp, decreasing = TRUE)[seq_len(minscreen)]] <- TRUE
  }

  names(whichVariable) <- colnames(X_df)

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
#' @examples
#' data("metabric", package = "SuperSurv")
#' dat <- metabric[1:40, ]
#' x_cols <- grep("^x", names(dat))[1:6]
#' X <- dat[, x_cols, drop = FALSE]
#'
#' screen.var(
#'   time = dat$duration,
#'   event = dat$event,
#'   X = X,
#'   keep_fraction = 0.5,
#'   minscreen = 2
#' )
#' @export
screen.var <- function(time, event, X, obsWeights = NULL,
                       keep_fraction = 0.5, minscreen = 2, ...) {

  X_mat <- .make_design_matrix(X)

  vars <- apply(X_mat, 2, function(col) stats::var(col, na.rm = TRUE))

  cutoff <- stats::quantile(vars, probs = 1 - keep_fraction, na.rm = TRUE)
  whichVariable <- (vars >= cutoff)
  whichVariable[is.na(whichVariable)] <- FALSE

  if (sum(whichVariable) < minscreen) {
    whichVariable <- rep(FALSE, length(vars))
    whichVariable[order(vars, decreasing = TRUE)[seq_len(minscreen)]] <- TRUE
  }
  names(whichVariable) <- colnames(X_mat)

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
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data("metabric", package = "SuperSurv")
#'   dat <- metabric[1:40, ]
#'   x_cols <- grep("^x", names(dat))[1:5]
#'   X <- dat[, x_cols, drop = FALSE]
#'
#'   screen.elasticnet(
#'     time = dat$duration,
#'     event = dat$event,
#'     X = X,
#'     alpha = 0.5,
#'     minscreen = 2,
#'     nfolds = 3,
#'     nlambda = 20
#'   )
#' }
#' @export
screen.elasticnet <- function(time, event, X, obsWeights = NULL, alpha = 0.5, minscreen = 2, nfolds = 10, nlambda = 100, ...) {

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
