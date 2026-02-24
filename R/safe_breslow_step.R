#' Calculate Baseline Survival using Breslow Estimator
#'
#' @param time Numeric vector of observed follow-up times.
#' @param event Numeric vector of event indicators (0=censored, 1=event).
#' @param risk_score Numeric vector of risk scores.
#' @param new.times Numeric vector of times at which to evaluate the baseline hazard.
#' 
#' @return A numeric vector representing the baseline survival probabilities 
#'   evaluated at the specified \code{new.times}.
#' @export
#' @keywords internal

safe_breslow_step <- function(time, event, risk_score, new.times) {
  
  # 1. Sort data by time
  ord <- order(time, -event)
  time_sorted <- time[ord]
  event_sorted <- event[ord]
  
  # Safety: Clamp huge scores to avoid Inf.
  safe_score <- pmin(risk_score[ord], 700) 
  risk_sorted <- exp(safe_score)
  
  # 2. Compute Risk Set Sums
  risk_pool <- rev(cumsum(rev(risk_sorted)))
  
  # 3. Compute Hazard contributions
  is_event <- event_sorted == 1
  if (!any(is_event)) return(rep(0, length(new.times)))
  
  event_t <- time_sorted[is_event]
  haz_contribution <- 1 / risk_pool[is_event]
  
  # 4. Aggregate ties
  unique_t <- unique(event_t)
  if (length(unique_t) < length(event_t)) {
    haz_contribution <- tapply(haz_contribution, event_t, sum)
    event_t <- as.numeric(names(haz_contribution))
  }
  
  # 5. Cumulative Hazard
  cum_haz <- cumsum(haz_contribution)
  
  # 6. Interpolate (Step Function)
  stats::approx(
    x = event_t, 
    y = cum_haz, 
    xout = new.times, 
    method = "constant", 
    yleft = 0, 
    yright = max(cum_haz),
    f = 0, 
    rule = 2,
    ties = mean
  )$y
}