#' Control parameters for Cross-Validation in SuperSurv
#'
#' @param V Number of folds. Default is 10.
#' @param stratifyCV Logical. If TRUE, ensures event rates are balanced across folds.
#' @param shuffle Logical. If TRUE, shuffles rows before splitting.
#' @param validRows Optional custom list of indices for folds.
#'
#' @return A list of CV parameters.
#' @export
SuperSurv.CV.control <- function(V = 10L, 
                                 stratifyCV = TRUE, 
                                 shuffle = TRUE, 
                                 validRows = NULL) {
  
  V <- as.integer(V)
  
  # Logic check: Can't have more folds than observations
  # (This will be double-checked in the main function once N is known)
  
  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("validRows must be a list of length V.")
    }
    if (!identical(V, length(validRows))) {
      # If user provided 5 folds of indices but set V=10, we must stop.
      V <- length(validRows)
      warning(paste("V adjusted to match length of validRows:", V))
    }
  }
  
  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}