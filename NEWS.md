## 2016-01-27

  * Handle `constant_item = FALSE` correctly in `wrangle()`

## 2016-01-24

  * Bugfix: group order could be scrambled in `ZZ` and `XX`

## 2016-01-18

  * Bugfix: adjusted trial and success counts (`n_vec`, `s_vec`), after calculation, could be associated with the wrong
    group names
  * Switch to numeric time variable in dgirt output

## 2016-01-10

  * Handle item variables (defensively) as follows:
    * A numeric item variable with two unique values or an ordered factor item variable with two observed levels
      represents a binary choice in which the higher value or level is "success" and the lower value or level "failure"
    * Generally, numeric and ordered-factor item variables represent ascending ordinal choices
    * Other classes of item variable (e.g. factor, character) result in an error
    * Print messages explaining the handling of each item variable
  * Remove dependency on `mcgv` in `plot_means()`

## 2016-01-06

  * Bugfix: mean item outcomes could be calculated incorrectly, inflating success counts toward trial counts.

## 2015-12-30

Version bump to 0.0.10.

  * Functionality:
    * Apply more descriptive names to `dgirt()` results using the variable names originally passed to `wrangle()` and
      the levels of factors.
    * `poststratify()` is safer and more flexible. It takes new arguments `strata`, `groups`, and `check_proportions`;
      see the documentation.
    * Specify the algorithm for CmdStan to use. `dgirt()` passes new argument `optimize_algorithm` to CmdStan if `method
      = "optimize"`; one of `"lbfgs"` (the default), `"bfgs"` and `"newton"`.
  * Documentation:
    * Switch to a README.Rmd that includes the "Getting Started" vignette content and drop the vignette.
  * Functions renamed:
    * `run_dgirt()` -> `dgirt()`
    * `format_dgirt()` -> `wrangle()`
  * New and renamed datasets:
    * `rstan_output`: example of `dgirt()` output for `method = "rstan"`
    * `optimize_output`: example of `dgirt()` output for `method = "optimize"`
    * `states` -> `state_opinion`
    * `state_targets` and `targets` -> `state_demographics`
  * Under the hood:
    * Switch to [assertthat](https://github.com/hadley/assertthat) package from ad-hoc `stop()` calls
    * Speed up `wrangle()`
    * Bugfixes
