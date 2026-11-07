#' Configure the optional CMG backend for two-way OLS
#'
#' Requires the experimental cmgr package. This backend accepts exactly two
#' categorical intercept fixed effects. It leaves native defaults unchanged.
#' @param fixef_tol Relative weighted fixed-effect moment tolerance. This is not
#'   numerically equivalent to the native MAP iterate-change tolerance.
#' @param fixef_maxiter Maximum PCG iterations per residualized column.
#' @param schur Enable exact low-degree elimination in the larger fixed-effect
#'   block (normally workers). Enabled by default; FALSE retains the full CMG graph.
#' @param tol,maxiter Compatibility aliases for fixef_tol and fixef_maxiter.
#'   Do not supply an alias and its corresponding canonical argument together.
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
