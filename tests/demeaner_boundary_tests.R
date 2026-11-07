# Public input shapes and native control boundaries.
library(fixest)
setFixest_nthreads(1)
setFixest_notes(FALSE)
set.seed(42)
a <- rep(1:500,each=6)
b <- rep(sample(1:30,500,TRUE),each=6)
move <- sample(length(a),150)
b[move] <- sample(1:30,150,TRUE)
d <- data.frame(a,b,x=rnorm(length(a)),z=rnorm(length(a)),y=rnorm(length(a)))
d$endog <- d$z+rnorm(nrow(d))
d$group <- d$a>250

# Unlike feols/MapDemeaner, native demean() has always allowed any positive tol.
raw <- as.matrix(d[c("y","x")])
f <- list(d$a,d$b)
ordinary <- demean(raw,f,tol=1e-13)
trusted <- demean(raw,f,slope.flag=c(0L,0L),slope.vars=list(0),weights=1,
                   fixef.algo=demeaning_algo(),tol=1e-13,im_confident=TRUE)
stopifnot(identical(ordinary,trusted),all(is.finite(ordinary)))

# Explicit tolerance/iteration overrides on a native prepared environment work.
env <- feols(y~x|a+b,d,fixef.tol=1e-3,only.env=TRUE)
resumed <- feols(env=env,fixef.tol=1e-10,fixef.iter=20000)
direct <- feols(y~x|a+b,d,fixef.tol=1e-10,fixef.iter=20000)
stopifnot(identical(coef(resumed),coef(direct)),
          identical(resid(resumed),resid(direct)),
          identical(resumed$iterations,direct$iterations))
env <- suppressWarnings(feols(y~x|a+b,d,fixef.tol=1e-10,fixef.iter=1,only.env=TRUE))
resumed <- feols(env=env,fixef.iter=20000)
stopifnot(identical(resid(resumed),resid(direct)))

# Constructor controls and split samples inherit the resolved settings.
configured <- feols(y~x|a+b,d,demeaner=MapDemeaner(fixef_tol=1e-10,fixef_maxiter=20000))
stopifnot(identical(resid(configured),resid(direct)))
splits <- feols(y~x|a+b,d,split=~group,fixef.tol=1e-10,fixef.iter=20000)
for(i in 1:2){
  single <- feols(y~x|a+b,d,subset=d$group==(i==2),fixef.tol=1e-10,fixef.iter=20000)
  stopifnot(identical(resid(splits[[i]]),resid(single)),
            identical(splits[[i]]$iterations,single$iterations))
}

# Compare documented vector, matrix and list/data.frame shapes with dense QR.
small <- d[seq_len(180),]
reference <- qr.resid(qr(model.matrix(~factor(a)+factor(b),small)),
                      as.matrix(small[c("y","x")]))
controls <- list(map=MapDemeaner(fixef_tol=1e-10))
if(requireNamespace("withinr",quietly=TRUE)){
  controls$withinr <- LsmrDemeaner(fixef_atol=1e-10,fixef_btol=1e-10)
}

f <- list(small$a,small$b)
for(name in names(controls)){
  ctl <- controls[[name]]
  vector <- demean(small$x,f,demeaner=ctl)
  frame <- demean(small[c("y","x")],f,as.matrix=TRUE,demeaner=ctl)
  one_frame <- demean(small["x"],f,as.matrix=TRUE,demeaner=ctl)
  list_result <- demean(unclass(small[c("y","x")]),f,as.matrix=TRUE,demeaner=ctl)
  stopifnot(identical(dim(vector),c(nrow(small),1L)),
            max(abs(vector[,1]-reference[,2]))<1e-7,
            max(abs(frame-reference))<1e-7,
            max(abs(one_frame[,1]-reference[,2]))<1e-7,
            max(abs(list_result-reference))<1e-7)
  if(!inherits(ctl,"fixest_demeaner_map")){
    env <- feols(y~x|a+b,d,demeaner=ctl,only.env=TRUE)
    error <- tryCatch(feols(env=env,fixef.tol=1e-10),error=identity)
    stopifnot(inherits(error,"error"),grepl("its constructor",conditionMessage(error),fixed=TRUE))
  }
  cat("Input/control boundary passed:",name,"\n")
}
