# Focused backend contract checks. Numerical comparisons use an independent
# dense dummy-variable QR reference on a small weighted disconnected graph.
library(fixest)
setFixest_nthreads(1)
set.seed(20261107)
n <- 900L
d <- data.frame(a=rep(1:90, each=10), x=rnorm(n), z=rnorm(n), w=exp(rnorm(n, sd=.4)))
d$b <- sample(1:10,n,TRUE) + 10L*(d$a>45)
d$c <- rep(1:3, length.out=n)
d$y <- .7*d$x+rnorm(90)[d$a]+rnorm(20)[d$b]+rnorm(n)
d$endog <- .8*d$z+rnorm(n)
d$y2 <- d$y+.4*rnorm(n)
fails <- function(expr, pattern=NULL){
  result <- tryCatch(force(expr), error=identity)
  stopifnot(inherits(result,"error"))
  if(!is.null(pattern)) stopifnot(grepl(pattern,conditionMessage(result),fixed=TRUE))
}

# Choosing MAP explicitly must preserve the ordinary native path exactly.
native <- feols(y~x|a+b,d,weights=~w,demeaned=TRUE)
map <- feols(y~x|a+b,d,weights=~w,demeaned=TRUE,demeaner=MapDemeaner())
stopifnot(identical(coef(native),coef(map)), identical(resid(native),resid(map)),
          identical(native$X_demeaned,map$X_demeaned),
          identical(native$y_demeaned,map$y_demeaned),
          identical(native$iterations,map$iterations))
map_tight <- MapDemeaner(fixef_tol=1e-10,fixef_maxiter=20000L)
native_tight <- feols(y~x|a+b,d,weights=~w,fixef.tol=1e-10,fixef.iter=20000)
chosen <- feols(y~x|a+b,d,weights=~w,demeaner=map_tight)
stopifnot(identical(coef(native_tight),coef(chosen)),identical(resid(native_tight),resid(chosen)))
fails(feols(y~x|a+b,d,demeaner=MapDemeaner(),fixef.tol=1e-8),"do not also provide")
fails(demean(d$x,list(d$a,d$b),demeaner=MapDemeaner(),iter=10),"do not also provide")
fails(feols(y~x|a+b,d,demeaner=structure(list(),class="fixest_demeaner")),"created by")
# Settable defaults follow the same interface and never alter existing defaults.
local({
  old <- options("fixest_estimation"); on.exit(options(old))
  poisson <- data.frame(y=c(1,3,2,4,1,2,3,5),x=c(0,1,2,3,0,1,2,3),a=rep(1:2,each=4))
  baseline <- fepois(y~x|a,poisson)
  setFixest_estimation(reset=TRUE,demeaner=map_tight)
  z <- feols(y~x|a+b,d,weights=~w)
  stopifnot(identical(coef(z),coef(chosen)),identical(resid(z),resid(chosen)))
  fails(setFixest_estimation(demeaner=MapDemeaner(),fixef.iter=20))
  implicit <- fepois(y~x|a,poisson)
  stopifnot(identical(coef(baseline),coef(implicit)),identical(resid(baseline),resid(implicit)))
})

controls <- list(map=map_tight)
if(requireNamespace("cmgr",quietly=TRUE)){
  controls$cmg <- CmgDemeaner(fixef_tol=1e-10,fixef_maxiter=10000L,schur=FALSE)
  controls$schur <- CmgDemeaner(fixef_tol=1e-10,fixef_maxiter=10000L,schur=TRUE)
  stopifnot(identical(controls$cmg,CmgDemeaner(tol=1e-10,maxiter=10000L)))
  fails(CmgDemeaner(fixef_tol=1e-8,tol=1e-8))
  fails(CmgDemeaner(fixef_maxiter=100,maxiter=100))
  fails(feols(y~x|a,d,demeaner=controls$cmg),"exactly two")
  fails(feols(y~x|a+b+c,d,demeaner=controls$cmg),"exactly two")
  fails(feols(y~x|a[z]+b,d,demeaner=controls$cmg),"exactly two")
}
if(requireNamespace("withinr",quietly=TRUE)){
  controls$withinr <- LsmrDemeaner(fixef_atol=1e-10,fixef_btol=1e-10,
                                  fixef_maxiter=10000L,preconditioner="additive")
  fails(feols(y~x|a[z]+b,d,demeaner=controls$withinr),"varying slopes")
  one <- feols(y~x|a,d,demeaner=controls$withinr)
  ref <- feols(y~x|a,d,fixef.tol=1e-10)
  stopifnot(identical(coef(one),coef(ref)),identical(resid(one),resid(ref)))
  three <- feols(y~x|a+b+c,d,weights=~w,demeaner=controls$withinr)
  oracle <- lm.wfit(model.matrix(~factor(a)+factor(b)+factor(c)+x,d),d$y,d$w)
  stopifnot(max(abs(resid(three)-oracle$residuals))<1e-7)
  fails(feols(y~x|a+b,d,demeaner=LsmrDemeaner(fixef_atol=1e-14,
       fixef_btol=1e-14,fixef_maxiter=1L,preconditioner="off")),"did not converge")
}
oracle <- lm.wfit(model.matrix(~factor(a)+factor(b)+x,d),d$y,d$w)
projection <- lm.wfit(model.matrix(~factor(a)+factor(b),d),as.matrix(d[c("y","x")]),d$w)$residuals
for(name in names(controls)){
  ctl <- controls[[name]]
  fit <- feols(y~x|a+b,d,weights=~w,demeaned=TRUE,vcov=~a+b,demeaner=ctl)
  stopifnot(max(abs(resid(fit)-oracle$residuals))<1e-7,
            abs(coef(fit)[["x"]]-tail(oracle$coefficients,1))<1e-8)
  dm <- demean(d[c("y","x")],list(d$a,d$b),weights=d$w,demeaner=ctl)
  stopifnot(max(abs(as.matrix(dm)-projection))<1e-7)
  matrix_fit <- feols.fit(d$y,as.matrix(d["x"]),d[c("a","b")],weights=d$w,demeaner=ctl)
  stopifnot(max(abs(resid(matrix_fit)-resid(fit)))<1e-7)
  saved <- feols(y~x|a+b,d,weights=~w,demeaned=TRUE)
  stopifnot(max(abs(demean(saved,demeaner=ctl)-projection))<1e-7)
  prepared <- feols(y~x|a+b,d,weights=~w,demeaner=ctl,only.env=TRUE)
  stopifnot(max(abs(coef(feols(env=prepared))-coef(fit)))<1e-8)
  fails(feols(env=prepared,demeaner=ctl),"creating that environment")
  fails(feglm(y~x|a+b,d,family="gaussian",demeaner=ctl),"OLS only")
  for(fml in list(y~x|a+b|endog~z, c(y,y2)~sw0(x)|a+b)){
    a <- feols(fml,d,weights=~w,demeaner=map_tight)
    b <- feols(fml,d,weights=~w,demeaner=ctl)
    aa <- if(inherits(a,"fixest_multi")) as.list(a) else list(a)
    bb <- if(inherits(b,"fixest_multi")) as.list(b) else list(b)
    for(i in seq_along(aa)) stopifnot(identical(obs(aa[[i]]),obs(bb[[i]])),
        max(c(0,abs(coef(aa[[i]])-coef(bb[[i]]))))<1e-8,
        max(abs(resid(aa[[i]])-resid(bb[[i]])))<1e-7)
  }
  cat("Backend contract passed:",name,"\n")
}
