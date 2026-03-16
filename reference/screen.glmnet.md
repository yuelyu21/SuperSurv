# GLMNET (Lasso) Screening

GLMNET (Lasso) Screening

## Usage

``` r
screen.glmnet(
  time,
  event,
  X,
  obsWeights = NULL,
  alpha = 1,
  minscreen = 2,
  nfolds = 10,
  nlambda = 100,
  ...
)
```

## Arguments

- time:

  Observed follow-up time.

- event:

  Observed event indicator.

- X:

  Training covariate data.frame.

- obsWeights:

  Observation weights.

- alpha:

  Penalty exponent (1 = lasso).

- minscreen:

  Minimum number of covariates to return. Defaults to 2.

- nfolds:

  Number of CV folds.

- nlambda:

  Number of penalty parameters.

- ...:

  Additional ignored arguments.

## Value

A logical vector of the same length as the number of columns in `X`,
indicating which variables passed the screening algorithm (`TRUE` to
keep, `FALSE` to drop).

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
  data("metabric", package = "SuperSurv")
  dat <- metabric[1:40, ]
  x_cols <- grep("^x", names(dat))[1:5]
  X <- dat[, x_cols, drop = FALSE]

  screen.glmnet(
    time = dat$duration,
    event = dat$event,
    X = X,
    alpha = 1,
    minscreen = 2,
    nfolds = 3,
    nlambda = 20
  )
}
#>    x0    x1    x2    x3    x4 
#>  TRUE  TRUE FALSE FALSE FALSE 
```
