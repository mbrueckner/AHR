library(testthat)
library(AHR)
	
test_that("Aalen-Johansen estimator reduces to Kaplan-Meier estimator for two-state model without recovery", {
    test.data <- function(n) {
        T <- rexp(n, 0.1)
        C <- runif(n, 0, 10)
        X <- pmin(T, C)
        D <- as.numeric(T <= C)
        V <- runif(n, 0, X/4) 

        status <- T <= C
        D[D == 0] <- "cens"
        
        data.frame(time=X, from=0, to=D, id=1:n, status=status)
    }
    
    data <- test.data(100)
    tra <- matrix(FALSE, nrow=2, ncol=2)
    tra[1, 2] <- TRUE
              
    times <- seq(0, 5, length.out=10) ##c(1.5, 3)

    data <- data[order(data$time),]

    data$to[1] <- 1
    data$status[1] <- TRUE
    
    ##S <- exp(-times/10)

    fs <- survfit(Surv(time, status) ~ 1, data=data)

    ## estimate cumulative incidence function for event type 1
    fit <- aj(sort(data$time), data, list(target="0 1", states=c("0", "1"), transitions=tra,
                                          censoring="cens", s=0, t="last", covariance=TRUE))

    
    f <- approxfun(fs$time, fs$surv, method="constant", yleft=1, rule=2, f=0)
    f2 <- approxfun(fs$time, fit$S, method="constant", yleft=1, rule=2, f=0)
    g <- approxfun(fs$time, (fs$surv * fs$std.err)^2, method="constant", yleft=0, rule=2, f=0)
    g2 <- approxfun(fs$time, fit$V, method="constant", yleft=0, rule=2, f=0)

    expect_true(all.equal(f(times), f2(times)))
    expect_true(all.equal(g(times), g2(times)))
})


test.etm <- function(j=1) {
    data <- data.frame(id=1:10, time=1:10, from=0, to=1, status=TRUE)

    tra <- matrix(FALSE, nrow=2, ncol=2)
    tra[1, 2] <- TRUE

    data$to[j] <- "cens"
    data$status[j] <- FALSE

    fit.km <- survfit(Surv(time, status) ~ 1, data=data)
    fit.etm <- etm(data, c("0","1"), tra, "cens", s=0, t="last", covariance=FALSE)

    data.frame(time=fit.km$time[data$status],
               km=fit.km$surv[data$status],         
               time2=as.numeric(names(fit.etm$est[1,2,data$status])),
               etm=1-fit.etm$est[1,2,data$status])
}
