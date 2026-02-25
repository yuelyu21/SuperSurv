# Predict method for SuperSurv fits

Obtains predicted survival probabilities from a fitted SuperSurv
ensemble.

## Usage

``` r
# S3 method for class 'SuperSurv'
predict(object, newdata, new.times, onlySL = FALSE, threshold = 1e-04, ...)
```

## Arguments

- object:

  A fitted object of class `SuperSurv`.

- newdata:

  A data.frame of new covariate values.

- new.times:

  A numeric vector of times at which to predict survival.

- onlySL:

  Logical. If TRUE, only uses models with weights \> threshold.

- threshold:

  Numeric. The weight threshold for onlySL.

- ...:

  Additional ignored arguments.

## Value

A list containing:

- `event.SL.predict`: A numeric matrix (rows = observations, columns =
  times) of the final predicted survival probabilities from the
  ensemble.

- `event.library.predict`: A 3D numeric array (observations x times x
  models) containing the individual survival predictions from each base
  learner.

- `cens.SL.predict`: A numeric matrix of the predicted censoring
  probabilities.

- `cens.library.predict`: A 3D numeric array of the individual censoring
  predictions.
