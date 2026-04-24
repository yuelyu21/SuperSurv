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

  Logical; if `TRUE`, compute perturbation-based standard errors,
  confidence intervals, and a Wald-type p-value. Defaults to `FALSE`.

- B:

  Integer giving the number of perturbation replicates when
  `inference = TRUE`. Defaults to `200`.

- seed:

  Optional integer seed for reproducibility of the perturbation
  procedure.

- ci_level:

  Numeric scalar in `(0, 1)` specifying the confidence level for the
  Wald-type confidence interval. Defaults to `0.95`.

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

- `inference`: Logical indicator for whether perturbation-based
  inference was requested.

- `B`: Number of perturbation replicates used when `inference = TRUE`;
  otherwise `NULL`.

- `SE_RMST`: Perturbation-based standard error of the RMST contrast;
  otherwise `NULL`.

- `CI_RMST`: Wald-type confidence interval for the RMST contrast;
  otherwise `NULL`.

- `z_value`: Wald-type test statistic; otherwise `NULL`.

- `p_value`: Two-sided Wald-type p-value; otherwise `NULL`.

- `perturb_reps`: Vector of perturbation replicate estimates; otherwise
  `NULL`.

## Details

For a binary grouping variable `trt_col`, the function predicts
counterfactual survival curves under `A = 1` and `A = 0` for every
individual in the supplied dataset, integrates each curve up to the
restriction time `tau`, and averages the resulting individual-level RMST
differences. The resulting contrast is generally interpreted as an
adjusted marginal contrast. When `trt_col` corresponds to a manipulable
intervention and additional identification assumptions hold, the same
standardized procedure may also support a causal interpretation.

If `inference = TRUE`, the function additionally performs a
perturbation-based inference procedure conditional on the fitted
`SuperSurv` model. In this implementation, the fitted learner library,
hyperparameters, base learners, and ensemble weights are held fixed, and
random positive weights are applied to the individual-level RMST
contrasts to estimate a perturbation-based standard error, Wald-type
confidence interval, and p-value.

The function uses the empirical distribution of the observed covariates
in `data` as the standardization distribution. RMST is evaluated
numerically from the predicted survival matrix using a left Riemann sum
over the supplied grid `times`.

The perturbation-based inference implemented here is conditional on the
fitted `SuperSurv` model. It does not re-tune hyperparameters, reselect
the learner library, or refit the base learners under each perturbation.
Instead, it perturbs the aggregation of the individual-level predicted
RMST contrasts. This yields a lightweight uncertainty quantification
procedure for the standardized RMST contrast given the final fitted
ensemble.

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
  tau = 100,
  inference = TRUE,
  B = 200,
  seed = 123
)

rmst_res$ATE_RMST
rmst_res$SE_RMST
rmst_res$CI_RMST
format.pval(rmst_res$p_value, digits = 3, eps = 1e-16)
} # }
```
