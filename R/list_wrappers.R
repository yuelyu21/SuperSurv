#' List Available Wrappers and Screeners in SuperSurv
#'
#' This function prints all built-in prediction algorithms and feature screening 
#' algorithms available in the \code{SuperSurv} package.
#'
#' @param what Character string. If \code{"both"} (default), lists both prediction 
#'   and screening functions. If \code{"surv"}, lists only prediction models. 
#'   If \code{"screen"}, lists only screening algorithms. Otherwise, lists all exports.
#' 
#' @return An invisible character vector containing the requested function names.
#' @export
list_wrappers <- function(what = "both") {
  
  # Grab all exported functions from the package
  everything <- sort(getNamespaceExports("SuperSurv"))
  
  # Identify models and screeners using regex
  models <- everything[grepl("^surv\\.", everything)]
  screens <- everything[grepl("^screen\\.", everything)]
  
  if (what == "both") {
    message("--- Prediction Models (surv.*) ---")
    print(models)
    message("\n--- Screening Algorithms (screen.*) ---")
    print(screens)
    invisible(c(models, screens))
    
  } else if (what == "surv") {
    message("--- Prediction Models (surv.*) ---")
    print(models)
    invisible(models)
    
  } else if (what == "screen") {
    message("--- Screening Algorithms (screen.*) ---")
    print(screens)
    invisible(screens)
    
  } else {
    message("--- All Exported Functions in SuperSurv ---")
    print(everything)
    invisible(everything)
  }
}