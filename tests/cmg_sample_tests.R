library(fixest)
if(requireNamespace("cmgr", quietly=TRUE)) local({
  setFixest_notes(FALSE)
  set.seed(20261109)
  a <- rep(1:120, each=12)
  b <- unlist(lapply(1:120,function(i)
    rep(c(i %% 6L + 1L, (i + 1L) %% 6L + 1L), 6) + 6L * (i > 60)))
  d <- data.frame(a, b, x=rnorm(length(a)), w=exp(rnorm(length(a), sd=.4)))
  worker_effect <- rnorm(120)[a]
  firm_effect <- rnorm(12)[b]
  d$y <- .8*d$x + worker_effect + firm_effect + rnorm(length(a))

  # FE-only columns should be excluded even when their original scale is large.
  # Tight native projection and independent weighted dummy QR are references.
  for(kind in c("worker", "firm", "additive")) for(scale in c(1, 1e5)){
    d$q <- scale * switch(kind, worker=worker_effect, firm=firm_effect,
                          additive=worker_effect+firm_effect)
    qr <- lm.wfit(model.matrix(~factor(a)+factor(b)+x+q,d),d$y,w=d$w)
    native <- feols(y~x+q|a+b,d,weights=~w,fixef.tol=1e-10,nthreads=1)
    stopifnot(is.na(qr$coefficients[["q"]]), "q" %in% native$collin.var)
    for(schur in c(FALSE,TRUE)){
      fit <- feols(y~x+q|a+b,d,weights=~w,nthreads=1,
                   demeaner=CmgDemeaner(schur=schur))
      stopifnot(identical(fit$collin.var,native$collin.var),
                abs(coef(fit)[["x"]]-qr$coefficients[["x"]]) < 1e-8,
                max(abs(resid(fit)-qr$residuals)) < 1e-7)
    }
  }

  # Cleaning the sample must precede low-degree elimination and preserve rows.
  d$q <- NULL
  d$y2 <- d$y+rnorm(nrow(d))
  d$y[2] <- NA_real_
  d$x[3] <- NA_real_
  d$w[4] <- 0
  d$a[nrow(d)] <- 121L
  d$b[nrow(d)] <- 13L
  native <- feols(y~x|a+b,d,weights=~w,fixef.tol=1e-10,nthreads=1)
  check_qr <- function(fit, response){
    sample <- d[obs(fit),]
    qr <- lm.wfit(model.matrix(~factor(a)+factor(b)+x,sample),
                   sample[[response]],w=sample$w)
    stopifnot(abs(coef(fit)[["x"]]-qr$coefficients[["x"]]) < 1e-8,
              max(abs(resid(fit)-qr$residuals)) < 1e-7)
  }
  for(schur in c(FALSE,TRUE)) for(threads in c(1L,2L)){
    fit <- feols(y~x|a+b,d,weights=~w,nthreads=threads,
                 demeaner=CmgDemeaner(fixef_tol=1e-10,schur=schur))
    stopifnot(identical(obs(fit),obs(native)),
              !any(c(2L,3L,4L,nrow(d)) %in% obs(fit)),
              max(abs(resid(fit)-resid(native))) < 1e-7)
    check_qr(fit,"y")
    # Confirm these cleaned samples actually take the elimination path.
    selected <- d[obs(fit),]
    projected <- cmgr::demean_fixest(selected$y,as.matrix(selected["x"]),
      fit$fixef_id,as.integer(fit$fixef_sizes),weights=selected$w,
      tol=1e-10,schur=schur,nthreads=threads)
    if(schur) stopifnot(projected$reduction[["eliminated_levels"]] > 0)
  }

  # Distinct outcome missingness removes a worker from one model and a firm
  # from the other. Each model must use its own compact fixed-effect IDs.
  d$y[d$a == 120L] <- NA_real_
  d$y2[d$b == 1L] <- NA_real_
  for(schur in c(FALSE,TRUE)){
    ctl <- CmgDemeaner(fixef_tol=1e-10,schur=schur)
    multiple <- as.list(feols(c(y,y2)~x|a+b,d,weights=~w,nthreads=2,demeaner=ctl))
    singles <- list(feols(y~x|a+b,d,weights=~w,nthreads=1,demeaner=ctl),
                    feols(y2~x|a+b,d,weights=~w,nthreads=1,demeaner=ctl))
    stopifnot(!identical(obs(multiple[[1]]),obs(multiple[[2]])),
              !120L %in% d$a[obs(multiple[[1]])],
              !1L %in% d$b[obs(multiple[[2]])])
    for(i in seq_along(singles)){
      stopifnot(identical(obs(multiple[[i]]),obs(singles[[i]])),
                max(abs(coef(multiple[[i]])-coef(singles[[i]]))) < 1e-9,
                max(abs(resid(multiple[[i]])-resid(singles[[i]]))) < 1e-9)
      check_qr(multiple[[i]],c("y","y2")[[i]])
    }
  }
  cat("CMG FE-collinearity, cleaned low-degree samples, multiple-outcome IDs and threads passed\n")
})
