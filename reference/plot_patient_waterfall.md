# Waterfall Plot for an Individual Patient

Waterfall Plot for an Individual Patient

## Usage

``` r
plot_patient_waterfall(shap_values, patient_index = 1, top_n = 10)
```

## Arguments

- shap_values:

  The output from explain_shap.\*

- patient_index:

  The row index of the patient to explain

- top_n:

  Number of features to show (default 10)

## Value

A `ggplot` object visualizing the SHAP values.
