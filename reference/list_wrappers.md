# List Available Wrappers and Screeners in SuperSurv

This function prints all built-in prediction algorithms and feature
screening algorithms available in the `SuperSurv` package.

## Usage

``` r
list_wrappers(what = "both")
```

## Arguments

- what:

  Character string. If `"both"` (default), lists both prediction and
  screening functions. If `"surv"`, lists only prediction models. If
  `"screen"`, lists only screening algorithms. Otherwise, lists all
  exports.

## Value

An invisible character vector containing the requested function names.

## Examples

``` r
list_wrappers()
#> --- Prediction Models (surv.*) ---
#>  [1] "surv.aorsf"       "surv.bart"        "surv.coxboost"    "surv.coxph"      
#>  [5] "surv.exponential" "surv.gam"         "surv.gbm"         "surv.glmnet"     
#>  [9] "surv.km"          "surv.loglogistic" "surv.lognormal"   "surv.parametric" 
#> [13] "surv.ranger"      "surv.rfsrc"       "surv.ridge"       "surv.rpart"      
#> [17] "surv.svm"         "surv.weibull"     "surv.xgboost"    
#> 
#> --- Screening Algorithms (screen.*) ---
#> [1] "screen.all"        "screen.elasticnet" "screen.glmnet"    
#> [4] "screen.marg"       "screen.rfsrc"      "screen.var"       
```
