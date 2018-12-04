library(testthat)
library(AHR)
     
test_check("AHR")

test_that("User supplied survival function estimator is the same as built-in KM estimator", {
T <- c(rexp(100, 1), rexp(100, 2))
C <- c(rexp(100, 1), rexp(100, 2))
time <- pmin(T, C)
status <- T <= C
trt <- rep(c(0,1), c(100, 100)) # treatment indicator

sfit <- function(times, data, param) {
    fit <- survfit(Surv(Y, D) ~ 1, data=data)
    f <- approxfun(fit$time, fit$surv, method="constant", f=0, yleft=1, rule=2)
    fv <- approxfun(fit$time, fit$std.err^2, method="constant", f=0, yleft=0, rule=2)

    S <- f(times)
    V <- fv(times)
        
    dlogS <- 1/S
    dlogS[S == 0] <- 0
    logV <- V * nrow(data)
    logV[S == 0] <- 0
    
    list(times=times, S=S, V=V, logV=logV)
}

fit1 <- ahrUser(2, Surv(time, status) ~ trt, data.frame(time=time, status=status, trt=trt), user.survfit=sfit, user.param=list())

fit2 <- ahrKM(2, Surv(time, status) ~ trt, data.frame(time=time, status=status, trt=trt), cov=TRUE)

expect_equal(fit1$surv.fit[[1]]$S, fit2$surv.fit[[1]]$S)
expect_equal(fit1$surv.fit[[2]]$S, fit2$surv.fit[[2]]$S)
expect_equal(fit1$surv.fit[[1]]$logV, fit2$surv.fit[[1]]$logV)
expect_equal(fit1$surv.fit[[2]]$logV, fit2$surv.fit[[2]]$logV)
})
