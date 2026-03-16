# Time-Dependent AUC and Integrated AUC

Evaluates the cumulative/dynamic time-dependent AUC and integrated AUC
(iAUC) using inverse probability of censoring weighting (IPCW).

## Usage

``` r
eval_timeROC(time, event, S_mat, times)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- S_mat:

  A numeric matrix of predicted survival probabilities.

- times:

  Numeric vector of evaluation times matching the columns of `S_mat`.

## Value

A list containing the `AUC_curve` at each time point, the `times`, and
the integrated AUC `iAUC`.

## Examples

``` r
 data("metabric", package = "SuperSurv")
 dat <- metabric[1:40, ]
 x_cols <- grep("^x", names(dat))[1:3]
 X <- dat[, x_cols, drop = FALSE]
 newX <- X[1:10, , drop = FALSE]
 times <- seq(50, 150, by = 50)

 fit <- surv.coxph(
   time = dat$duration,
   event = dat$event,
   X = X,
   newdata = newX,
   new.times = times,
   obsWeights = rep(1, nrow(dat)),
   id = NULL
 )

 eval_timeROC(
   time = dat$duration[1:10],
   event = dat$event[1:10],
   S_mat = fit$pred,
   times = times
 )
#> $AUC_curve
#> [1] 1.0000000 0.5000000 0.7192982
#> 
#> $times
#> [1]  50 100 150
#> 
#> $iAUC
#> [1] 0.6798246
#> 
```
