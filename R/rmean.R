#' Estimate difference of restricted mean survival based on (weighted) Kaplan-Meier estimates
#' of the survival functions in each group.
#'
#' @title rmeanDiff
#' @param L time-limit specifying up to which time restricted mean will be calculated
#' @param formula an object of class '"formula"' specifying the conditional survival model
#' @param data data frame containing the variables in formula
#' @return An object of class '"rmd"', i.e. a list containing:
#'  \item{L}{time limit, i.e. restricted mean up to time L is calculated}
#'  \item{rmean1}{restricted mean in group 1}
#'  \item{rmean2}{restricted mean in group 2}
#'  \item{rmean.diff}{estimated restricted mean difference}
#'  \item{var.rmean1}{an estimate of the asymptotic variance of the restricted mean in group 1}
#'  \item{var.rmean2}{an estimate of the asymptotic variance of the restricted mean in group 2}
#'  \item{var.rmean.diff}{an estimate of the asymptotic variance of the restricted mean difference}
#'  \item{Z.rmean}{the standardized test statistic for testing rmean.diff=0}
#'  \item{p.value}{p-value corresponding to Z.rmean}
#' @export
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100))
#' fit <- rmeanDiff(2, formula=Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z))
rmeanDiff <- function(L, formula, data) {        
    if(!is.null(formula)) data <- parseFormula(formula, data)
    
    grps <- levels(data$Trt)
    if(nlevels(data$Trt) != 2) stop("Difference in restricted means can only be estimated in two-sample setting")

    times <- getTimes(L, data)

    data1 <- data[data$Trt == grps[1], ]
    data2 <- data[data$Trt == grps[2], ]
    
    n <- c(nrow(data1), nrow(data2))

    param <- list(alpha=1, start=0, var=TRUE, cov=FALSE, left.limit=FALSE)

    fit1 <- wkm(times, data1, param)
    fit2 <- wkm(times, data2, param)

    rmeanDiffFit(L, times, fit1, fit2, n, n/sum(n), TRUE)    
}
                      
#' Estimate difference of restricted mean survival (based on ahr object as returned by ahr)
#'
#' This function is usefull if the function 'ahr' has already been called, since the survival estimates
#' in the object returned by 'ahr' can be reused.
#'
#' @title rmeanDiff.ahr
#' @param ahr.obj object of class '"ahr"'
#' @return An object of class '"rmd"', i.e. a list containing:
#'  \item{L}{time limit, i.e. restricted mean up to time L is calculated}
#'  \item{rmean1}{restricted mean in group 1}
#'  \item{rmean2}{restricted mean in group 2}
#'  \item{rmean.diff}{estimated restricted mean difference}
#'  \item{var.rmean1}{an estimate of the asymptotic variance of the restricted mean in group 1}
#'  \item{var.rmean2}{an estimate of the asymptotic variance of the restricted mean in group 2}
#'  \item{var.rmean.diff}{an estimate of the asymptotic variance of the restricted mean difference}
#'  \item{Z.rmean}{the standardized test statistic for testing rmean.diff=0}
#'  \item{p.value}{p-value corresponding to Z.rmean}
#' @export
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100))
#' fit <- avgHR(2, data.frame(Y=Y, D=D, Z=Z), formula=Surv(Y, D) ~ Z)
#' rmd <- rmeanDiff.ahr(fit)
rmeanDiff.ahr <- function(ahr.obj) {
    if(class(ahr.obj) != "ahr") stop("Input must be of class 'ahr'")
    if(nlevels(ahr.obj$groups) != 2) stop("Difference in restricted means can only be estimated in two-sample setting")

    rmeanDiffFit(ahr.obj$L, ahr.obj$times, ahr.obj$surv.fit[[1]], ahr.obj$surv.fit[[2]], ahr.obj$n, ahr.obj$p, ahr.obj$log.iis)
}

rmeanDiffFit <- function(L, times, fit1, fit2, n, p, log.iis) {

    snt <- 1:length(times)
    
    ## restricted means
    rmean1 <- stepIntegrate(fit1$S, times)
    rmean2 <- stepIntegrate(fit2$S, times)
    rmean.diff <- rmean1 - rmean2
    
    if((!is.null(fit1$COV) && !is.null(fit2$COV)) || log.iis) {
        if(log.iis) {
            v1 <- fit1$V / p[1]
            v2 <- fit2$V / p[2]
            
            ##tmp <- simplify2array(lapply(snt, function(l) stepIntegrate(v1[pmin(snt, l)] + v2[pmin(snt, l)], times)))
            tmp1 <- simplify2array(lapply(snt, function(l) stepIntegrate(v1[pmin(snt, l)], times)))
            tmp2 <- simplify2array(lapply(snt, function(l) stepIntegrate(v2[pmin(snt, l)], times)))
        } else {
              rho1 <- fit1$COV / p[1]
              rho2 <- fit2$COV / p[2]
              
              ##tmp <- simplify2array(lapply(snt, function(l) stepIntegrate(rho1[,l] + rho2[,l], times)))
              tmp1 <- simplify2array(lapply(snt, function(l) stepIntegrate(rho1[,l], times)))
              tmp2 <- simplify2array(lapply(snt, function(l) stepIntegrate(rho2[,l], times)))
          }
        
        ##var.rmean.diff <- stepIntegrate(tmp, times) / sqrt(sum(n))
        var.rmean1 <- stepIntegrate(tmp1, times) / sqrt(sum(n))
        var.rmean2 <- stepIntegrate(tmp2, times) / sqrt(sum(n))
        var.rmean.diff <- var.rmean1 + var.rmean2
        
        ## standardized difference
        Z <- rmean.diff / sqrt(var.rmean.diff)
          
        res <- list(rmean1=rmean1, rmean2=rmean2, rmean.diff=rmean.diff, var.rmean.diff=var.rmean.diff, var.rmean1=var.rmean1, var.rmean2=var.rmean2, Z.rmean=Z, p.value=1-pnorm(Z), L=L)
    } else res <- list(rmean.diff=rmean.diff, L=L)

    class(res) <- "rmd"
    res
}

#' Print rmd object
#'
#' @title print.rmd
#' @param x an object of class '"rmd"'.
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#' @method print rmd
#' @export
print.rmd <- function(x, digits=3, ...) {
    cat("\n")
    cat("\t Restricted mean difference", paste0("(L=", x$L, ")"), "\n")
    cat("\n\n")

    cat(paste0("Restricted mean in group 1: ", round(x$rmean1, digits), " (Std.dev: ", round(sqrt(x$var.rmean1), digits), ")\n"))
    cat(paste0("Restricted mean in group 2: ", round(x$rmean2, digits), " (Std.dev: ", round(sqrt(x$var.rmean2), digits), ")"))

    cat("\n\n")

    cat(paste0("Restricted mean difference: ", round(x$rmean.diff, digits), " (Std.dev: ", round(sqrt(x$var.rmean.diff), digits), ", Z: ", round(x$Z.rmean, digits), ", p: ", round(x$p.value, digits), ")"))

    cat("\n")
    
}
