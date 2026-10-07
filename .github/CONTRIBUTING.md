# Release documentation

`NEWS.md` is the versioned changelog published by pkgdown at
<https://yuelyu21.github.io/SuperSurv/news/index.html>.

For every package release:

1. Add a new `# SuperSurv x.y.z` entry at the top of `NEWS.md`. Retain earlier
   release entries and identify the previous public release being compared.
2. Include a `## Backward-Compatibility Summary` section. Describe changes to
   arguments, defaults, return values, numerical results, dependencies, and
   deprecations, including any action users need to take. When no such changes
   exist, state that explicitly.
3. Summarize new features and fixes separately. Verify each entry against the
   release source; do not infer unchanged numerical results from passing checks.
4. Confirm that `DESCRIPTION` matches the intended source release and that the
   changelog covers it. Keep private review correspondence out of the website.
5. Run `pkgdown::init_site()` followed by `pkgdown::build_news()` to preview the
   changelog, and review its content.
6. After the release changes are merged into the default branch, verify that the
   pkgdown workflow succeeds and that the public Changelog displays the new entry.

The existing pkgdown workflow rebuilds and publishes the website on default-branch
updates. Release notes must be written with each release; the workflow publishes
them but does not infer compatibility changes automatically.
