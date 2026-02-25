# Time-Dependent AUC and Integrated AUC

Evaluates the cumulative/dynamic Time-Dependent AUC and Integrated AUC
(iAUC) using the `timeROC` package with IPCW adjustment.

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

A list containing the `AUC` at each time point and the `iAUC`.
