# Estimate an Adjusted Marginal RMST Contrast

Computes a covariate-adjusted marginal contrast on the restricted mean
survival time (RMST) scale using standardization (g-computation) based
on a fitted `SuperSurv` model.

## Usage

``` r
estimate_marginal_rmst(
  fit,
  data,
  trt_col,
  times,
  tau,
  inference = FALSE,
  B = 200,
  seed = NULL,
  ci_level = 0.95
)
```

## Arguments

- fit:

  A fitted object of class `"SuperSurv"`.

- data:

  A `data.frame` containing the covariates used for standardization,
  including the binary grouping variable specified by `trt_col`.

- trt_col:

  Character string giving the name of the binary grouping variable in
  `data`. The variable is set to 1 and 0, respectively, to generate the
  two standardized prediction regimes.

- times:

  Numeric vector of prediction time points corresponding to the
  evaluation grid used for survival prediction.

- tau:

  Numeric scalar giving the restriction horizon for RMST. Must not
  exceed `max(times)`.

- inference:

  Deprecated compatibility argument. Formal inference is not provided;
  it requires a validated procedure accounting for model-estimation
  uncertainty. Must remain `FALSE`.

- B, seed, ci_level:

  Deprecated compatibility arguments; ignored when `inference = FALSE`.

## Value

A list containing:

- `ATE_RMST`: The estimated adjusted marginal RMST contrast
  \\\widehat{\Delta}\_{RMST}(\tau)\\.

- `mean_RMST_Treated`: The average predicted RMST under `A = 1`.

- `mean_RMST_Control`: The average predicted RMST under `A = 0`.

- `tau`: The restriction horizon used for integration.

- `patient_rmst_treated`: Vector of individual-level predicted RMST
  values under `A = 1`.

- `patient_rmst_control`: Vector of individual-level predicted RMST
  values under `A = 0`.

- `patient_delta_rmst`: Vector of individual-level predicted RMST
  contrasts.

- `inference`: Always `FALSE`; retained for compatibility.

## Details

For a binary grouping variable `trt_col`, the function predicts
counterfactual survival curves under `A = 1` and `A = 0` for every
individual in the supplied dataset, integrates each curve up to the
restriction time `tau`, and averages the resulting individual-level RMST
differences. The resulting contrast is generally interpreted as an
adjusted marginal contrast. When `trt_col` corresponds to a manipulable
intervention and additional identification assumptions hold, the same
standardized procedure may also support a causal interpretation.

The function uses the empirical distribution of the observed covariates
in `data` as the standardization distribution. RMST is evaluated
numerically from the predicted survival matrix using a left Riemann sum
over the supplied grid `times`.

The returned contrast is a model-based point estimate. Formal
uncertainty quantification would need to account for nuisance
estimation, tuning, cross-validation, learner fitting, and
ensemble-weight estimation, for example through a validated
full-pipeline refitting procedure.

## Examples

``` r
if (FALSE) { # \dontrun{
data("metabric", package = "SuperSurv")
x_cols <- grep("^x", names(metabric), value = TRUE)
X <- metabric[, x_cols]
new.times <- seq(10, 150, by = 10)

fit <- SuperSurv(
  time = metabric$duration,
  event = metabric$event,
  X = X,
  newdata = X,
  new.times = new.times,
  event.library = c("surv.coxph", "surv.rfsrc"),
  cens.library = c("surv.coxph"),
  control = list(saveFitLibrary = TRUE),
  nFolds = 3
)

rmst_res <- estimate_marginal_rmst(
  fit = fit,
  data = metabric,
  trt_col = "x4",
  times = new.times,
  tau = 100
)

rmst_res$ATE_RMST
} # }
```
