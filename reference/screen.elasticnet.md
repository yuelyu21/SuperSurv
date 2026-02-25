# Elastic Net Screening Algorithm

This screening algorithm uses
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html) to
select covariates. Unlike LASSO (`alpha = 1`), which drops correlated
features, Elastic Net (`alpha = 0.5` by default) shrinks correlated
groups of features together, making it ideal for selecting entire
biological pathways.

## Usage

``` r
screen.elasticnet(
  time,
  event,
  X,
  obsWeights,
  alpha = 0.5,
  minscreen = 2,
  nfolds = 10,
  nlambda = 100,
  ...
)
```

## Arguments

- time:

  Numeric vector of observed follow-up times.

- event:

  Numeric vector of event indicators (1 = event, 0 = censored).

- X:

  Training covariate data.frame or matrix.

- obsWeights:

  Numeric vector of observation weights.

- alpha:

  Numeric penalty exponent for `glmnet`. Defaults to 0.5 (Elastic Net).

- minscreen:

  Integer. Minimum number of covariates to return. Defaults to 2.

- nfolds:

  Integer. Number of folds for cross-validation. Defaults to 10.

- nlambda:

  Integer. Number of penalty parameters to search over. Defaults to 100.

- ...:

  Additional arguments passed to `screen.glmnet`.

## Value

A logical vector of the same length as the number of columns in `X`,
indicating which variables passed the screening algorithm (`TRUE` to
keep, `FALSE` to drop).
