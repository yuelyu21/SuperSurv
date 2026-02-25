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
