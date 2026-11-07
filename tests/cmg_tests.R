library(fixest)
if(requireNamespace("cmgr",quietly=TRUE)) for(schur in c(FALSE, TRUE)) local({
setFixest_nthreads(1)
set.seed(97012)
n <- 2400L
d <- data.frame(a=rep(1:120,each=20),x=rnorm(n),x2=rnorm(n),z=rnorm(n),w=exp(rnorm(n,sd=.6)))
d$b <- sample(1:15,n,TRUE)+15L*(d$a>60)
d$y <- .7*d$x-.3*d$x2+rnorm(120)[d$a]+rnorm(30)[d$b]+rnorm(n)
d$y2 <- d$y+.3*rnorm(n);d$endog <- .8*d$z+.3*d$x+rnorm(n)
ctl <- CmgDemeaner(tol=1e-10, schur=schur)
check <- function(fml,...) {
    cat("Checking",deparse(fml),"arguments",names(list(...)),"\n")
    a <- feols(fml,d,...,fixef.tol=1e-10)
    b <- feols(fml,d,...,demeaner=ctl)
    aa <- if(inherits(a,"fixest_multi")) as.list(a) else list(a)
    bb <- if(inherits(b,"fixest_multi")) as.list(b) else list(b)
    for(i in seq_along(aa)) {
        x<-aa[[i]];y<-bb[[i]]
        stopifnot(identical(obs(x),obs(y)),identical(names(coef(x)),names(coef(y))),
                  max(abs(coef(x)-coef(y)))<1e-8,
                  max(abs(resid(x)-resid(y)))<1e-7,
                  max(abs(vcov(x)-vcov(y)))<1e-8,
                  max(abs(predict(x)-predict(y)),na.rm=TRUE)<1e-7)
    }
    invisible(b)
}
check(y~x+x2|a+b,vcov=~a+b)
fit<-check(y~x+x2|a+b,weights=~w,vcov=~a+b)
check(y~x+x2|a+b|endog~z,weights=~w)
check(c(y,y2)~x+sw0(x2)|a+b,weights=~w)
check(y~x+x2|a+b,weights=~w,offset=~z)
check(y~x+x2|a+b,weights=~w,split=~I(a>30))
# Direct calls preserve formula expressions for prediction on new data.
o1 <- feols(y~x+x2|a+b,d,weights=~w,offset=~z,fixef.tol=1e-10)
o2 <- feols(y~x+x2|a+b,d,weights=~w,offset=~z,demeaner=ctl)
stopifnot(max(abs(predict(o1,newdata=d)-predict(o2,newdata=d)))<1e-7)
f1 <- feols(y~1|a+b,d,weights=~w,fixef.tol=1e-10)
f2 <- feols(y~1|a+b,d,weights=~w,demeaner=ctl)
stopifnot(max(abs(resid(f1)-resid(f2)))<1e-7)
m1 <- feols.fit(d$y,as.matrix(d[c("x","x2")]),d[c("a","b")],weights=d$w,fixef.tol=1e-10)
m2 <- feols.fit(d$y,as.matrix(d[c("x","x2")]),d[c("a","b")],weights=d$w,demeaner=ctl)
stopifnot(max(abs(coef(m1)-coef(m2)))<1e-8,max(abs(resid(m1)-resid(m2)))<1e-7)
oracle<-lm.wfit(model.matrix(~factor(a)+factor(b)+x+x2,d),d$y,w=d$w)
stopifnot(max(abs(tail(oracle$coefficients,2)-coef(fit)))<1e-8)
saved <- feols(y~x+x2|a+b,d,weights=~w,demeaned=TRUE)
stopifnot(max(abs(demean(saved,demeaner=ctl)-demean(saved)))<1e-5)
prepared <- feols(y~x+x2|a+b,d,weights=~w,demeaner=ctl,only.env=TRUE)
stopifnot(max(abs(coef(feols(env=prepared))-coef(fit)))<1e-8)
stopifnot(inherits(try(feols(env=prepared,demeaner=ctl),silent=TRUE),"try-error"))
d$collinear<-d$x
fit<-check(y~x+x2+collinear|a+b,weights=~w)
stopifnot("collinear" %in% fit$collin.var)
d$y[2]<-NA;d$x[3]<-NA;d$w[4]<-0;d$a[n]<-121;d$b[n]<-31
check(y~x+x2|a+b,weights=~w)
dd<-na.omit(d[d$w>0,])
dm<-demean(dd[c("y","x")],list(dd$a,dd$b),weights=dd$w,demeaner=ctl)
ref<-lm.wfit(model.matrix(~factor(a)+factor(b),dd),as.matrix(dd[c("y","x")]),w=dd$w)$residuals
stopifnot(max(abs(as.matrix(dm)-ref))<1e-8)
fails<-function(expr)stopifnot(inherits(try(expr,silent=TRUE),"try-error"))
bad <- ctl; bad$schur <- NA
fails(feols(y~x|a+b,d,demeaner=bad))
# Old saved configurations remain compatible and use the full graph.
old <- CmgDemeaner(); old$schur <- NULL
old_fit <- feols(y~x+x2|a+b,d,weights=~w,demeaner=old)
off_fit <- feols(y~x+x2|a+b,d,weights=~w,demeaner=CmgDemeaner(schur=FALSE))
stopifnot(identical(coef(old_fit),coef(off_fit)),identical(resid(old_fit),resid(off_fit)))
fails(feols(y~x,d,demeaner=ctl))
fails(feols(y~x|a,d,demeaner=ctl))
fails(feols(y~x|a+b+I(a>30),d,demeaner=ctl))
fails(feols(y~x|a[z]+b,d,demeaner=ctl))
fails(feols(y~x|a+b,d,demeaner=ctl,fixef.tol=1e-8))
fails(feglm(y~x|a+b,d,family="gaussian",demeaner=ctl))
cat("fixest CMG integration tests passed; schur=",schur,"\n")
})

# Exercise actual low-degree elimination in complete fits. Compare to the
# unchanged CMG-only path; dense dummy QR independently checks weighted OLS.
if(requireNamespace("cmgr",quietly=TRUE)) local({
set.seed(20261108)
a <- rep(1:120,each=6)
b <- unlist(lapply(1:120,function(i) {
  path <- if(i %% 5L == 0L) 1:6 else rep(c(i %% 6L+1L,(i+1L) %% 6L+1L),3)
  path + 6L*(i>60)
}))
d <- data.frame(a,b,x=rnorm(length(a)),z=rnorm(length(a)),w=exp(rnorm(length(a))))
d$y <- .8*d$x+rnorm(120)[a]+rnorm(12)[b]+rnorm(length(a))
d$endog <- d$z+rnorm(length(a)); d$y2 <- d$y+rnorm(length(a))
compare <- function(fml,...) {
  off <- feols(fml,d,...,demeaned=TRUE,demeaner=CmgDemeaner(tol=1e-10,schur=FALSE))
  on <- feols(fml,d,...,demeaned=TRUE,demeaner=CmgDemeaner(tol=1e-10,schur=TRUE))
  aa <- if(inherits(off,"fixest_multi")) as.list(off) else list(off)
  bb <- if(inherits(on,"fixest_multi")) as.list(on) else list(on)
  for(i in seq_along(aa)) stopifnot(identical(obs(aa[[i]]),obs(bb[[i]])),
    max(abs(coef(aa[[i]])-coef(bb[[i]])))<1e-9,
    max(abs(resid(aa[[i]])-resid(bb[[i]])))<1e-9,
    max(abs(vcov(aa[[i]])-vcov(bb[[i]])))<1e-9)
  invisible(on)
}
f <- compare(y~x|a+b,weights=~w,vcov=~a+b)
qr <- lm.wfit(model.matrix(~factor(a)+factor(b)+x,d),d$y,w=d$w)
stopifnot(max(abs(resid(f)-qr$residuals))<1e-9,
          abs(coef(f)[['x']]-tail(qr$coefficients,1))<1e-9)
compare(y~x|b+a,weights=~w)
compare(y~x|a+b|endog~z,weights=~w)
compare(c(y,y2)~x|a+b,weights=~w,split=~I(a>60))
compare(y~x|a+b,weights=~w,offset=~z)
cat("Actual Schur elimination in weighted OLS/IV/multiple/split/offset fits passed\n")
})
