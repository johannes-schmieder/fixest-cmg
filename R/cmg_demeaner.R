#' Configure the optional two-way CMG backend
#'
#' Configures explicit CMG residualization for two categorical intercept fixed
#' effects in OLS or [demean()]. Requires the optional experimental cmgr package.
#' The native backend remains the default.
#'
#' Supply the result as `demeaner = CmgDemeaner()` to [feols()], [feols.fit()] or
#' [demean()]. Selecting it per call is recommended in this preview: a global
#' selector set with [setFixest_estimation()] applies to every `feols` call, so
#' calls without exactly two categorical intercept effects fail by design.
#'
#' Put convergence controls in this constructor; do not also supply `fixef.tol`,
#' `fixef.iter` or `fixef.algo`. With `only.env = TRUE`, set the backend and its
#' controls when creating the environment; supplying the backend or native
#' controls when reusing that environment is an error.
#'
#' Other numbers of fixed effects, varying slopes, stepwise fixed-effect
#' specifications and GLM requests fail with an error. Solver failure also raises
#' an error, with no fallback to another solver. An interrupt raises an ordinary R
#' error, which `try()` and error handlers can catch.
#'
#' IV diagnostics and named fixed-effect extraction use native demeaning
#' internally, even on a CMG fit. Their residualization does not necessarily reach
#' the accuracy of the CMG projection.
#'
#' **Exact worker elimination.** A worker with few employers does not need to be
#' part of the graph that CMG solves. Given the firm effects, that worker's effect
#' is the worker's weighted mean of the outcome net of those firm effects.
#' Substituting this into the system removes the worker exactly and connects the
#' worker's firms directly. After the smaller system is solved, the worker's
#' effect is recovered.
#'
#' This is applied by default, as preprocessing inside cmgr. It works on the
#' fixed-effect block with more levels (normally workers; the first block on a
#' tie) and eliminates only levels with at most five distinct neighbors. Repeated
#' worker-firm observations are combined before counting. Higher-degree levels
#' stay in the solved graph. Nothing is sampled or approximated, no observations
#' are dropped, and CMG's options and PCG implementation are unchanged.
#'
#' `schur = FALSE` is an optional comparison and troubleshooting switch. It
#' disables this algebraic reduction, so CMG solves the full worker-firm graph for
#' the same model and observations. It does not select another solver or estimator.
#'
#' **Tolerance.** `fixef_tol` bounds a normalized weighted fixed-effect moment
#' recomputed from the returned observations. It differs from the native stopping
#' tolerance, so equal values do not mean equal accuracy, and it is not a guarantee
#' of coefficient accuracy on ill-conditioned graphs.
#'
#' @param fixef_tol Numeric scalar, finite and strictly between zero and one.
#'   Relative weighted fixed-effect moment tolerance; see Details.
#' @param fixef_maxiter Positive whole number. Maximum PCG iterations per
#'   residualized column.
#' @param schur `TRUE` (default) or `FALSE`. Whether to apply exact worker
#'   elimination before the CMG solve; see Details.
#' @param tol,maxiter Compatibility aliases for `fixef_tol` and `fixef_maxiter`.
#'   Supplying an alias together with its canonical argument is an error.
#' @return A configuration of class `c("fixest_demeaner_cmg", "fixest_demeaner")`
#'   with elements `fixef_tol`, `fixef_maxiter` and `schur`.
#' @seealso [MapDemeaner()], [LsmrDemeaner()]
#' @examples
#' if(requireNamespace("cmgr", quietly = TRUE)) {
#'   fit <- feols(mpg ~ wt | cyl + gear, mtcars, demeaner = CmgDemeaner())
#'   full_graph <- feols(mpg ~ wt | cyl + gear, mtcars, demeaner = CmgDemeaner(schur = FALSE))
#' }
#' @export
CmgDemeaner = function(fixef_tol = 1e-8, fixef_maxiter = 10000L, schur = TRUE,
                       tol = NULL, maxiter = NULL){
  if(!requireNamespace("cmgr", quietly = TRUE)){
    stop("CmgDemeaner() requires the optional 'cmgr' package.")
  }
  args = list(schur = schur)
  if(!missing(fixef_tol)) args$fixef_tol = fixef_tol
  if(!missing(fixef_maxiter)) args$fixef_maxiter = fixef_maxiter
  if(!missing(tol)) args$tol = tol
  if(!missing(maxiter)) args$maxiter = maxiter
  do.call(cmgr::CmgDemeaner, args)
}

#' @export
print.fixest_demeaner_cmg = function(x, ...){
  cat("<fixest_demeaner: CMG>\n")
  cat("  fixef_tol: ", if(is.null(x$fixef_tol)) x$tol else x$fixef_tol, "\n", sep = "")
  cat("  fixef_maxiter: ", if(is.null(x$fixef_maxiter)) x$maxiter else x$fixef_maxiter, "\n", sep = "")
  cat("  schur: ", if(is.null(x$schur)) TRUE else isTRUE(x$schur), "\n", sep = "")
  invisible(x)
}
