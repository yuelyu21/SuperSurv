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
#>  [1] "surv.aorsf"          "surv.bart"           "surv.coxboost"      
#>  [4] "surv.coxph"          "surv.coxtime"        "surv.deephit"       
#>  [7] "surv.deepsurv"       "surv.exponential"    "surv.flexsurvreg"   
#> [10] "surv.flexsurvspline" "surv.gam"            "surv.gbm"           
#> [13] "surv.glmnet"         "surv.grf"            "surv.km"            
#> [16] "surv.loglogistic"    "surv.lognormal"      "surv.mboost"        
#> [19] "surv.parametric"     "surv.ranger"         "surv.rfsrc"         
#> [22] "surv.ridge"          "surv.rpart"          "surv.survPen"       
#> [25] "surv.svm"            "surv.weibull"        "surv.xgboost"       
#> 
#> --- Screening Algorithms (screen.*) ---
#> [1] "screen.all"        "screen.elasticnet" "screen.glmnet"    
#> [4] "screen.marg"       "screen.rfsrc"      "screen.var"       
```
