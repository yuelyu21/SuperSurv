## JSS editor priorities

This repository is being revised for resubmission to the *Journal of Statistical Software (JSS)*.

The revision should directly address these package-facing editor comments:

- strengthen use of R classes and methods,
- provide and maintain `print.SuperSurv()` and `summary.SuperSurv()`,
- preserve `predict.SuperSurv()` and `coef.SuperSurv()`,
- provide clear user-facing accessors instead of encouraging direct `$` access to fitted-object internals,
- keep user-facing examples, vignettes, and pkgdown aligned with the revised object interface,
- document how new learners and screeners can be added,
- preserve statistical behavior unless a change is explicitly required for correctness.

## Current package interface expectations

The fitted `SuperSurv` object should behave like a mature R model object.

Expected user-facing interface includes:
- `print.SuperSurv()`
- `summary.SuperSurv()`
- `predict.SuperSurv()`
- `coef.SuperSurv()`
- `event_weights()`
- `censor_weights()`
- `learner_names()`
- `training_variables()`
- `selected_variables()`

Do not reintroduce user-facing examples that rely on undocumented direct list access such as `fit$event.coef`.

## Extensibility expectations

The package should present itself as an extensible framework, not just a fixed collection of wrappers.

Maintain clear documentation for:
- how learner wrappers are discovered,
- required wrapper inputs,
- required wrapper return structure,
- prediction matrix requirements,
- how direct survival-curve learners differ from risk-score learners,
- how screeners are defined.

If modifying wrapper-related docs, keep the custom-wrapper vignette practical and concise.

## Current known constraints

- The internal fitted-object structure may remain list-based.
- Avoid major refactors of internal storage unless clearly necessary.
- Do not destabilize prediction behavior just to simplify internals.
- Risk-score calibration helpers may remain internal if documented appropriately.
- If a new accessor would require object-structure redesign, prefer documenting the limitation unless the change is clearly worth it for JSS.

## Current remaining priorities

Highest remaining priorities in this repo:
1. keep docs/examples/vignettes synchronized with the redesigned object interface,
2. keep package/manuscript terminology consistent,
3. support final verification and version synchronization for resubmission.
