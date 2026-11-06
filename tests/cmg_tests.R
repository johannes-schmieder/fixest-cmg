library(fixest)
if(requireNamespace("cmgr",quietly=TRUE)) local({
setFixest_nthreads(1)
set.seed(97012)
n <- 2400L
d <- data.frame(a=rep(1:120,each=20),x=rnorm(n),x2=rnorm(n),z=rnorm(n),w=exp(rnorm(n,sd=.6)))
d$b <- sample(1:15,n,TRUE)+15L*(d$a>60)
d$y <- .7*d$x-.3*d$x2+rnorm(120)[d$a]+rnorm(30)[d$b]+rnorm(n)
d$y2 <- d$y+.3*rnorm(n);d$endog <- .8*d$z+.3*d$x+rnorm(n)
ctl <- CmgDemeaner(tol=1e-10)
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
fails(feols(y~x,d,demeaner=ctl))
fails(feols(y~x|a,d,demeaner=ctl))
fails(feols(y~x|a+b+I(a>30),d,demeaner=ctl))
fails(feols(y~x|a[z]+b,d,demeaner=ctl))
fails(feols(y~x|a+b,d,demeaner=ctl,fixef.tol=1e-8))
fails(feglm(y~x|a+b,d,family="gaussian",demeaner=ctl))
cat("fixest CMG integration tests passed\n")
})
