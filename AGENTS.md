# AGENTS.md

## Project
This repository contains the `SuperSurv` R package and related JSS revision materials for the manuscript:
**"SuperSurv: A Unified Framework for Machine Learning Ensembles in Survival Analysis"**

The immediate goal is a **high-quality redesign and revision** suitable for resubmission to the *Journal of Statistical Software (JSS)*.

## Primary goals
When working in this repository, prioritize the following:

1. Make the fitted `SuperSurv` object behave like a mature R model object.
2. Improve package extensibility and developer documentation.
3. Eliminate unnecessary direct `$` access to internals in examples, vignettes, and replication code.
4. Make manuscript replication complete, transparent, and easy to run.
5. Keep the manuscript, package version, examples, and replication materials synchronized.
6. Preserve statistical behavior unless a change is explicitly requested.

## Non-goals unless explicitly requested
Do not:
- redesign the statistical method itself unless needed for correctness,
- silently change numerical behavior,
- add major new methodology unrelated to the editor comments,
- remove supported functionality just to simplify code.

## Working style
Before making large edits:
1. Inspect the current code and summarize the relevant files.
2. Make a short implementation plan.
3. Prefer small, reviewable commits or diffs.

For substantial tasks, prefer **plan first, then edit**.

## Code standards
- Follow idiomatic R package structure.
- Prefer exported accessors and S3 methods over encouraging users to inspect internal lists.
- Keep user-facing APIs stable where possible.
- Use roxygen2 documentation for exported functions and methods.
- Update tests when behavior changes.
- Update examples/vignettes if APIs change.
- Avoid introducing unnecessary dependencies.
- Keep functions focused and readable.
- Preserve backward compatibility when reasonably possible.

## Required package improvements
The package should support, at minimum, the following user-facing interface:

- `print.SuperSurv()`
- `summary.SuperSurv()`
- `predict.SuperSurv()` (already exists or should be preserved)
- at least one accessor such as:
  - `coef.SuperSurv()` for ensemble weights, or
  - `weights.SuperSurv()`

If additional accessor functions are helpful, prefer clear names such as:
- `event_weights()`
- `censor_weights()`
- `learner_names()`
- `eval_times()`

## Object design principles
The fitted object may remain list-based internally, but users should not need to rely on undocumented internals.

Preferred approach:
- keep internal structure mostly intact if possible,
- add class-aware methods and accessors,
- ensure printing is concise,
- ensure summary output is informative.

## Replication principles
All outputs shown in the manuscript should be reproducible from the replication materials.

Preferred replication structure:
- one master script or one master R Markdown/Quarto file,
- clear README with how to run,
- figures/tables written to predictable output folders.

If manuscript examples currently depend on exploratory scripts, consolidate them.

## Manuscript revision principles
When editing the manuscript:
- preserve the core positioning of SuperSurv as a continuous-time survival Super Learner framework,
- clarify why this is not simply a thin extension of the `SuperLearner` package,
- explicitly discuss package extensibility,
- explain current scope choices such as right-censored data only,
- refer to R software as “packages,” not “libraries,” in polished prose.

## Specific editorial issues to address
Target the following JSS concerns explicitly:
- missing print/summary methods,
- insufficient use of classes and methods,
- direct `$` access in replication and examples,
- incomplete replication materials,
- unclear extensibility for new base learners,
- version mismatch between submission and CRAN,
- vignette ordering,
- comparison to broader software ecosystem and to `SuperLearner`,
- missing formal package references,
- absence of formula interface discussion,
- rationale for current wrapper subset,
- scope limitation to right-censoring.

## Testing and validation
After edits:
1. Run package checks if feasible.
2. Run affected examples/tests.
3. Check that new methods print sensible output.
4. Check that vignette/manuscript code no longer depends on undocumented internals where avoidable.
5. Flag any unresolved reproduction gaps clearly.

## Communication style
When reporting back:
- summarize what changed,
- list files modified,
- note any remaining open issues,
- identify anything that still requires human judgment.

## Repository map guidance
If both package and manuscript folders are present:
- treat the package folder as source of truth for APIs,
- update manuscript code examples only after package interface changes are settled,
- keep replication scripts aligned with the revised package API.

## Preferred workflow
For a new task:
1. inspect,
2. plan,
3. implement in small steps,
4. run checks,
5. summarize diff and any follow-up suggestions.
