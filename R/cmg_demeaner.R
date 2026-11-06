#' Configure the optional CMG backend for two-way OLS
#'
#' Requires the experimental cmgr package. This backend accepts exactly two
#' categorical intercept fixed effects. It leaves native defaults unchanged.
#' @param tol Relative weighted fixed-effect moment tolerance. This is not
#'   numerically equivalent to the native MAP iterate-change tolerance.
#' @param maxiter Maximum PCG iterations per residualized column.
#' @export
CmgDemeaner = function(tol = 1e-8, maxiter = 10000L){
  if(!requireNamespace("cmgr", quietly = TRUE)){
    stop("CmgDemeaner() requires the optional 'cmgr' package.")
  }
  cmgr::CmgDemeaner(tol, maxiter)
}

fixest_check_cmg_demeaner = function(demeaner){
  if(!inherits(demeaner, "fixest_demeaner_cmg")){
    stop("Argument 'demeaner' must be NULL or a CmgDemeaner() configuration.")
  }
  if(!requireNamespace("cmgr", quietly = TRUE)){
    stop("The CMG demeaner requires the optional 'cmgr' package.")
  }
  # Revalidate even when a configuration was modified or unserialized.
  cmgr::CmgDemeaner(demeaner$tol, demeaner$maxiter)
  invisible(NULL)
}

fixest_demean = function(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q,
                         fe_id_list, table_id_I, slope_flag_Q, slope_vars_list,
                         r_init, nthreads, algo_extraProj, algo_iter_warmup,
                         algo_iter_projAfterAcc, algo_iter_grandAcc,
                         demeaner = NULL){
  if(is.null(demeaner)){
    return(cpp_demean(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q,
                     fe_id_list, table_id_I, slope_flag_Q, slope_vars_list,
                     r_init, nthreads, algo_extraProj, algo_iter_warmup,
                     algo_iter_projAfterAcc, algo_iter_grandAcc))
  }
  fixest_check_cmg_demeaner(demeaner)
  if(length(r_nb_id_Q) != 2L || any(slope_flag_Q != 0L)){
    stop("CMG currently requires exactly two categorical intercept fixed effects.")
  }
  # fixest uses scalar one as its FE-only X sentinel; cmgr has an explicit
  # scalar-zero sentinel at its public low-level boundary.
  if(length(X_raw) == 1L) X_raw = 0
  if(!is.list(y) && length(y) == 1L) y = 0
  cmgr::demean_fixest(y, X_raw, fe_id_list, r_nb_id_Q,
                     weights = if(length(r_weights) == 1L) NULL else r_weights,
                     tol = demeaner$tol, maxiter = demeaner$maxiter,
                     nthreads = nthreads)
}
