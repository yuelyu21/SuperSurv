# SuperSurv 0.1.9

Changes since the previous CRAN release, 0.1.7.

## Backward-Compatibility Summary

The default survival transformation remains exponential. The following changes
can affect existing analyses or function calls.

### Fitting and prediction

* Direct simplex optimization can change ensemble weights relative to normalized
  nonnegative least squares (NNLS).
* Corrected tied-event calibration can change calibrated survival predictions and
  ensemble weights on data with tied event times.
* Right-continuous censoring endpoint corrections can change fitted weights at
  grid points coinciding with observed censoring times.
* Models can now be fitted without supplying prediction data or prediction times.
  Existing calls that supply both arguments remain supported.
* Product-limit conversion is opt-in; exponential conversion remains the default.
* `metalearner = "entropy"` remains temporarily available with a deprecation
  warning. Use `"brier"` or `"logloss"` for new analyses.
* Single-patient and single-time forest predictions retain their documented
  matrix dimensions.

### Evaluation and interpretation

* Corrected reverse Kaplan-Meier tie handling can change evaluation metrics on
  data with tied event and censoring times.
* RMST point contrasts remain available. The standard error, confidence interval,
  and p-value are no longer returned; `inference = TRUE` now produces an error.
* The optional SHAP backend is now `kernelshap`, and `explain_kernel()` requires
  an explicit `eval_time`. Explanations target the fitted ensemble's event
  probability at that time rather than a composite of native learner scores.
  The output table structure is retained, but earlier numerical SHAP results
  are not comparable. Recompute explanations with an explicitly chosen horizon.

## Documentation

* Shortened vignette builds for routine CRAN checks. The final vignette-runtime
  changes do not alter the fitting, prediction, or evaluation code.

