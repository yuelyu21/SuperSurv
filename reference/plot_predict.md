# Plot Predicted Survival Curves

Generates a ggplot2 step-function plot of the predicted survival
probabilities over time for one or more specific patients.

## Usage

``` r
plot_predict(preds, eval_times, patient_idx = 1)
```

## Arguments

- preds:

  A list containing predictions, specifically the object returned by the
  `predict.SuperSurv` function.

- eval_times:

  Numeric vector of times at which predictions were evaluated.

- patient_idx:

  Integer vector. The row indices of the patients in the test dataset
  whose survival curves you want to plot. Defaults to 1.

## Value

A `ggplot` object showing the survival curves.
