#' Demeaning backend constructors
#'
#' `MapDemeaner()` selects fixest's native method of alternating projections.
#' `LsmrDemeaner()` selects the optional `withinr` LSMR backend.
#'
#' @param fixef_tol Convergence tolerance for the MAP backend.
#' @param fixef_maxiter Maximum number of fixed-effect iterations.
#' @param extraProj,iter_warmup,iter_projAfterAcc,iter_grandAcc Arguments passed
#'   to [demeaning_algo()] for the MAP backend.
#' @param fixef_atol,fixef_btol Absolute tolerances for the LSMR backend. The
#'   current `withinr` API accepts one tolerance, so fixest passes the larger
#'   of the two values.
#' @param preconditioner LSMR preconditioner selection. Character values
#'   `"auto"`, `"additive"`, `"off"`, and `"diagonal"` are supported. Built
#'   `withinr` preconditioners and `withinr::AdditiveSchwarz()` configurations
#'   are passed through to `withinr`.
#' @param local_size Optional local reorthogonalization window size, or `NULL`.
#'
#' @return A `fixest_demeaner` object accepted by `feols()` and `demean()`.
#'
#' @export
MapDemeaner = function(fixef_tol = 1e-6, fixef_maxiter = 10000,
                       extraProj = 0, iter_warmup = 15,
                       iter_projAfterAcc = 40, iter_grandAcc = 4){

  check_arg(fixef_tol, "numeric scalar GT{(10000*.Machine$double.eps)}")
  check_arg(fixef_maxiter, "integer scalar GE{1}")
  check_arg("integer scalar", extraProj, iter_warmup, iter_projAfterAcc, iter_grandAcc)

  res = list(
    fixef_tol = fixef_tol,
    fixef_maxiter = as.integer(fixef_maxiter),
    fixef.algo = demeaning_algo(
      extraProj = extraProj,
      iter_warmup = iter_warmup,
      iter_projAfterAcc = iter_projAfterAcc,
      iter_grandAcc = iter_grandAcc,
      internal = TRUE
    )
  )

  class(res) = c("fixest_demeaner_map", "fixest_demeaner")
  res
}

#' @rdname MapDemeaner
#' @export
LsmrDemeaner = function(fixef_atol = 1e-6, fixef_btol = 1e-6,
                        fixef_maxiter = 10000, preconditioner = "auto",
                        local_size = NULL){

  if(!fixest_withinr_supports_lsmr_options()){
    stop("LsmrDemeaner() requires the optional 'withinr' package with ",
         "LsmrOptions() and within_solve_batch() or solve_batch().")
  }

  check_arg(fixef_atol, fixef_btol, "numeric scalar GT{0}")
  check_arg(fixef_maxiter, "integer scalar GE{1}")

  if(!is.null(local_size)){
    check_arg(local_size, "integer scalar GE{1}")
    local_size = as.integer(local_size)
  }

  preconditioner = fixest_check_lsmr_preconditioner(preconditioner)

  res = list(
    fixef_atol = fixef_atol,
    fixef_btol = fixef_btol,
    fixef_maxiter = as.integer(fixef_maxiter),
    preconditioner = preconditioner,
    local_size = local_size
  )

  class(res) = c("fixest_demeaner_lsmr", "fixest_demeaner")
  res
}

fixest_check_lsmr_preconditioner = function(preconditioner){

  if(is.null(preconditioner)){
    return(NULL)
  }

  if(is.character(preconditioner) && length(preconditioner) == 1L && !is.na(preconditioner)){
    return(preconditioner)
  }

  if(inherits(preconditioner, "within_preconditioner") ||
     inherits(preconditioner, "within_additive_schwarz")){
    return(preconditioner)
  }

  stop("Argument `preconditioner` must be a character scalar, NULL, ",
       "a withinr preconditioner, or a withinr AdditiveSchwarz configuration.")
}

fixest_check_demeaner = function(demeaner){

  if(is.null(demeaner)){
    return(invisible(NULL))
  }

  if(!inherits(demeaner, "fixest_demeaner") ||
     !any(vapply(c("fixest_demeaner_map", "fixest_demeaner_lsmr",
                   "fixest_demeaner_cmg"), function(cl) inherits(demeaner, cl), logical(1)))){
    stop("Argument `demeaner` must be NULL or created by MapDemeaner(), ",
         "LsmrDemeaner() or CmgDemeaner().")
  }

  invisible(NULL)
}

fixest_resolve_demeaner = function(demeaner, fixef.tol, fixef.iter, fixef.algo,
                                   explicit_demeaning_args = character()){

  if(!is.null(demeaner) && length(explicit_demeaning_args) > 0){
    pblm = paste0("`", explicit_demeaning_args, "`", collapse = ", ")
    stop("When `demeaner` is provided, do not also provide ",
         pblm, ". ",
         "Put these controls inside MapDemeaner(), LsmrDemeaner() or CmgDemeaner().")
  }

  if(is.null(fixef.algo)){
    fixef.algo = demeaning_algo(internal = TRUE)
  }

  if(is.null(demeaner)){
    demeaner = MapDemeaner(
      fixef_tol = fixef.tol,
      fixef_maxiter = fixef.iter,
      extraProj = fixef.algo$extraProj,
      iter_warmup = fixef.algo$iter_warmup,
      iter_projAfterAcc = fixef.algo$iter_projAfterAcc,
      iter_grandAcc = fixef.algo$iter_grandAcc
    )
  } else {
    fixest_check_demeaner(demeaner)
    # Revalidate modified or serialized configurations. Preserve old CMG
    # configurations from before the common fixef_* argument names.
    if(inherits(demeaner, "fixest_demeaner_cmg")){
      demeaner = CmgDemeaner(
        fixef_tol = if(is.null(demeaner$fixef_tol)) demeaner$tol else demeaner$fixef_tol,
        fixef_maxiter = if(is.null(demeaner$fixef_maxiter)) demeaner$maxiter else demeaner$fixef_maxiter,
        schur = if(is.null(demeaner$schur)) TRUE else demeaner$schur
      )
    } else if(inherits(demeaner, "fixest_demeaner_lsmr")){
      demeaner = do.call(LsmrDemeaner, unclass(demeaner))
    } else {
      algo = demeaner$fixef.algo
      demeaner = MapDemeaner(demeaner$fixef_tol, demeaner$fixef_maxiter,
                            algo$extraProj, algo$iter_warmup,
                            algo$iter_projAfterAcc, algo$iter_grandAcc)
    }
  }

  demeaner
}

fixest_demeaner_tol = function(demeaner){

  if(inherits(demeaner, "fixest_demeaner_lsmr")){
    return(max(demeaner$fixef_atol, demeaner$fixef_btol))
  }

  demeaner$fixef_tol
}

fixest_demeaner_iter = function(demeaner){
  demeaner$fixef_maxiter
}

fixest_demeaner_algo = function(demeaner){

  if(inherits(demeaner, "fixest_demeaner_map")){
    return(demeaner$fixef.algo)
  }

  demeaning_algo(internal = TRUE)
}

fixest_lsmr_preconditioner_arg = function(demeaner){

  preconditioner = demeaner$preconditioner

  if(is.character(preconditioner)){
    preconditioner_lower = tolower(preconditioner)

    if(preconditioner_lower == "auto"){
      return("additive")
    }

    if(preconditioner_lower %in% c("additive", "off", "diagonal")){
      return(preconditioner_lower)
    }

    warning("Unsupported LSMR preconditioner \"", preconditioner,
            "\"; using the withinr default.", call. = FALSE)
    return(NULL)
  }

  preconditioner
}

fixest_withinr_supports_lsmr_options = function(){
  requireNamespace("withinr", quietly = TRUE) &&
    "LsmrOptions" %in% getNamespaceExports("withinr") &&
    any(c("within_solve_batch", "solve_batch") %in% getNamespaceExports("withinr"))
}

#' @export
print.fixest_demeaner_map = function(x, ...){
  cat("<fixest_demeaner: MAP>\n")
  cat("  fixef_tol: ", x$fixef_tol, "\n", sep = "")
  cat("  fixef_maxiter: ", x$fixef_maxiter, "\n", sep = "")
  invisible(x)
}

#' @export
print.fixest_demeaner_lsmr = function(x, ...){
  preconditioner = x$preconditioner
  if(is.character(preconditioner)){
    preconditioner = paste0('"', preconditioner, '"')
  } else if(is.null(preconditioner)){
    preconditioner = "NULL"
  } else {
    preconditioner = class(preconditioner)[1]
  }

  cat("<fixest_demeaner: LSMR>\n")
  cat("  fixef_atol: ", x$fixef_atol, "\n", sep = "")
  cat("  fixef_btol: ", x$fixef_btol, "\n", sep = "")
  cat("  fixef_maxiter: ", x$fixef_maxiter, "\n", sep = "")
  cat("  preconditioner: ", preconditioner, "\n", sep = "")
  cat("  local_size: ", if(is.null(x$local_size)) "NULL" else x$local_size, "\n", sep = "")
  invisible(x)
}

# Explicit dispatcher following the interface proposed in fixest PR #647.
# fixest owns sample preparation; optional packages own numerical projection.
fixest_demean = function(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q,
                         fe_id_list, table_id_I, slope_flag_Q, slope_vars_list,
                         r_init, nthreads, fixef.algo, demeaner = NULL){
  fixest_check_demeaner(demeaner)
  if(is.null(demeaner) || inherits(demeaner, "fixest_demeaner_map") ||
     (inherits(demeaner, "fixest_demeaner_lsmr") && length(r_nb_id_Q) == 1L)){
    return(cpp_demean(y, X_raw, r_weights, iterMax, diffMax, r_nb_id_Q,
                     fe_id_list, table_id_I, slope_flag_Q, slope_vars_list,
                     r_init, nthreads,
                     algo_extraProj = fixef.algo$extraProj,
                     algo_iter_warmup = fixef.algo$iter_warmup,
                     algo_iter_projAfterAcc = fixef.algo$iter_projAfterAcc,
                     algo_iter_grandAcc = fixef.algo$iter_grandAcc))
  }

  # Native demeaning accepts lists and bare vectors as X. Optional backends
  # require an observation-by-column matrix; preserve scalar sentinels.
  if(is.list(X_raw)){
    X_raw = do.call(cbind, X_raw)
  } else if(is.numeric(X_raw) && is.null(dim(X_raw)) && length(X_raw) > 1L){
    X_raw = matrix(X_raw, ncol = 1L)
  }

  if(inherits(demeaner, "fixest_demeaner_cmg")){
    if(length(r_nb_id_Q) != 2L || any(slope_flag_Q != 0L)){
      stop("CMG currently requires exactly two categorical intercept fixed effects.")
    }
    # fixest uses scalar one as its FE-only X sentinel.
    if(length(X_raw) == 1L) X_raw = 0
    if(!is.list(y) && length(y) == 1L) y = 0
    return(cmgr::demean_fixest(y, X_raw, fe_id_list, r_nb_id_Q,
              weights = if(length(r_weights) == 1L) NULL else r_weights,
              tol = demeaner$fixef_tol, maxiter = demeaner$fixef_maxiter,
              nthreads = nthreads, schur = demeaner$schur))
  }

  fixest_demean_lsmr(y, X_raw, r_weights, fe_id_list, slope_flag_Q, demeaner)
}

fixest_demean_lsmr = function(y, X_raw, r_weights, fe_id_list, slope_flag_Q,
                              demeaner){
  if(any(slope_flag_Q != 0L)){
    stop("LsmrDemeaner() does not support varying slopes. Use MapDemeaner().")
  }
  is_y_list = is.list(y)
  has_y = is_y_list || (is.numeric(y) && length(y) > 1L)
  has_X = is.numeric(X_raw) && length(X_raw) > 1L
  y_mat = if(is_y_list) do.call(cbind, y) else if(has_y) as.matrix(y) else NULL
  X_mat = if(has_X) as.matrix(X_raw) else NULL
  n_y = NCOL(y_mat)
  n_X = NCOL(X_mat)
  Y = cbind(X_mat, y_mat)
  categories = do.call(cbind, fe_id_list)
  options = withinr::LsmrOptions(
    tol = max(demeaner$fixef_atol, demeaner$fixef_btol),
    maxiter = as.integer(demeaner$fixef_maxiter), local_size = demeaner$local_size)
  # within PR #84 renamed solve_batch to within_solve_batch. Support the
  # installable older R-binding preview as well as the new bindings.
  entry = if("within_solve_batch" %in% getNamespaceExports("withinr"))
    "within_solve_batch" else "solve_batch"
  result = getExportedValue("withinr", entry)(
    categories, Y, options = options,
    weights = if(length(r_weights) == 1L) NULL else r_weights,
    preconditioner = fixest_lsmr_preconditioner_arg(demeaner))
  if(!all(result$converged)){
    stop("withinr did not converge for all variables. Increase fixef_maxiter ",
         "or use MapDemeaner().")
  }
  dm = result$demeaned
  if(has_X){
    X_demean = dm[, seq_len(n_X), drop = FALSE]
    colnames(X_demean) = colnames(X_mat)
  } else {
    X_demean = 0
  }
  y_demean = if(!has_y) 0 else if(is_y_list)
    lapply(n_X + seq_len(n_y), function(i) dm[, i]) else dm[, n_X + 1L]
  list(X_demean = X_demean, y_demean = y_demean,
       iterations = as.integer(result$iterations), means = 0)
}
