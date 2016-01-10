#' Weighted Kaplan-Meier estimator with discrete time-independent covariate
#'
#' @title wkm
#' @param times a vector of evaluation times
#' @param data data frame containing the variables in formula
#' @param param list of parameters containing:
#' start: time of interim analysis (estimation of response rates is based only on data accruing after time 'start'), only if recruitment times are supplied (vector R in data) (default=0)
#' alpha: fractional parameter (default=1)
#' var: if TRUE (default) calculate variance estimate
#' cov: if FALSE (default) do not calculate covariance matrix estimate
#' left.limit: if TRUE calculate left-continuous estimates, else calculate right-continuous estimates
#' @param formula an object of class '"formula"' specifying the conditional survival model (only discrete covariates supported)
#' @return an object of class '"wkm"'
#' @details This function calculates the weighted Kaplan-Meier estimator for the survival function with weights based on a discrete time-independent covariate as described in Murray/Tsiatis (1996).
#' The survival probabilities are evaluated at each entry in the vector \code{times}. The data frame \code{data} must either contain the variable in \code{formula} or, if \code{formula} is \code{NULL},
#' the variables \code{V} (left-truncation time), \code{Y} (censored failure time), \code{D} (censoring indicator), \code{W} (stratification variable) and optionally \code{R} (recruitment time).
#' If \code{var} is \code{TRUE} then an estimate of the asmyptotic variance is calculated for each entry in vector \code{times}. If \code{cov} is \code{TRUE} then the \code{n x n} asymptotic
#' covariance matrix is estimated, where \code{n} is the length of vector \code{times}. If \code{left.limit} is \code{TRUE} then a left-continuous estimate of the survival function is calculated instead
#' of a right-continuous estimate (default). If recruitment times are supplied in \code{data} as a variable named \code{R} then the weights are estimated using observations from patients recruited after
#' calendar time \code{start} only.
#' @references S.~Murray and A.~A. Tsiatis. Nonparametric survival estimation using prognostic longitudinal covariates. \emph{Biometrics}, 52(1):137--151, Mar. 1996.
#' @export
wkm <- function(times, data, param=list(start=0, alpha=1, var=TRUE, cov=FALSE, left.limit=FALSE), formula=NULL) {

    if(is.null(param)) param <- list(start=0, alpha=1, var=TRUE, cov=FALSE, left.limit=FALSE)
    else {
        start <- param$start
        alpha <- param$alpha    
        var <- param$var
        cov  <- param$cov   
        left.limit <- param$left.limit
    }
    
    if(!is.null(formula)) data <- parseFormula(formula, data, one.sample=TRUE)
    ## if is.null(formula) assume that the variables in data are named V,Y,D,W

    strata <- levels(factor(data$W))
    n.strata <- length(strata)

    n.times <- length(times)
    n <- length(data$Y)
    
    S <- numeric(n.times)  
    P <- numeric(n.strata)
    n.atrisk <- numeric(n.times)

    if(var) {
        logV <- 0
        V <- 0
    } else {
        logV <- NULL
        V <- NULL
    }

    if(n.times == 1) cov <- FALSE
    
    if(cov) COV <- matrix(0, nrow=n.times, ncol=n.times)
    else {
        COV <- NULL
        logCOV <- NULL
    }

    ## second stage data (if start=0, then data2 == data)
    if(!is.null(data$R)) data2 <- data[data$R >= start, ]
    else data2 <- data

    ## set left-truncation times to 0 if not supplied by user
    if(is.null(data$V)) data$V <- rep.int(0, n)
    
    for(s in 1:n.strata) {
        s.data <- data[data$W == strata[s], ]

        if(nrow(s.data) == 0) stop("Empty stratum!")

        fit <- fastkm(s.data$Y, s.data$D, s.data$V, left.limit, times)

        ## total number-at-risk
        n.atrisk <- n.atrisk + fit$n.atrisk
        
        ## unbiased response rates
        p <- length(data2$W[data2$W == strata[s]]) / nrow(data2)
        stopifnot(!is.nan(p))
        
        P[s] <- p
        
        fs <- fit$surv
        tfs <- p * fs

        S <- S + tfs
        
        ## fit$variance (Greenwood) is variance of -log(\hat{S}), but we need variance of -\sqrt{ns}log(\hat{S})
        fvar <- fit$variance * nrow(s.data)
        fvar[fit$surv == 0] <- 0

        if(cov) { ## calculate covariance matrix
            ## exploit independent increments structure! (cov(X(s), X(t)) = var(X(min(s,t))))
            v <- matrix(fvar, nrow=n.times, ncol=n.times, byrow=FALSE) 
            v[lower.tri(v)] <- t(v)[lower.tri(v)]

            dfs <- diag(fs)
            
            COV <- COV + p * dfs %*% v %*% dfs + tfs %*% t(fs)            
        }
        
        if(var) { ## calculate variance
            V <- V + tfs * fs * (1 + fvar)
        }
    }

    Sa.log <- 1/(S^alpha)
    Sa.log[S == 0] <- 0
    
    if(cov) {
        COV <- COV - S %*% t(S)
        
        if(alpha != 1) {
            Sa <- S^(alpha - 1)
            Sa[S == 0] <- 0  ## if alpha - 1 < 0 then S == 0 implies V = NA but should be V = 0
            ss <- diag(Sa)
            COV <- alpha^2 * ss %*% COV %*% ss
        }
        
        ssl <- diag(Sa.log)
        logCOV <- ssl %*% COV %*% ssl

        ##V <- diag(COV)
        ##logV <- diag(logCOV)        
    }
    
    if(var) { 
        V <- V - S^2
        
        if(alpha != 1) {
            Sa <- S^(alpha - 1)
            Sa[S == 0] <- 0  ## if alpha - 1 < 0 then S == 0 implies V = NA but should be V = 0
            V <- alpha^2 * Sa^2 * V
        }
        
        logV <- Sa.log^2 * V
    }

    obj <- list(times=times, alpha=alpha, start=start, p=P, S=S^alpha, COV=COV, logCOV=logCOV, logV=logV, V=V, n.atrisk=n.atrisk)
    class(obj) <- "wkm"
    obj
}


################################################
### WKM with time-dependent discrete covariate
################################################
wkmSelect <- function(path, times, data) {
    if(is.null(path)) data
    else {
        col <- length(path)
        sel <- (data$Y > times[col]) & (data[, col+2] == path[col])
        data[sel, ]
    }
}

## depth-first traversal of all possible covariate paths
wkmTraverse <- function(path, k, tt, times, data) {
  s.data <- wkmSelect(path, times, data)
  n <- dim(s.data)[1]
  if(is.null(n)) c(0, 0) ## no subjects on this covariate path
  else if(tt <= times[length(path)+1]) { ## leaf
    c(n, n * fastkm(s.data$Y, s.data$D, s.data$V, FALSE, tt)$surv[1])
  } else {
    x <- sapply(0:k, function(i) wkmTraverse(c(path, i), k, tt, times, s.data))
    tmp <- sum(x[2,]) / sum(x[1,])
    if(is.null(path)) tmp  ## root
    else c(n, n * fastkm(s.data$Y, s.data$D, s.data$V, FALSE, times[length(path)])$surv[1] * tmp)
  }
}

wkmTime <- function(times, k, data) {
  sapply(times, function(tt) wkmTraverse(NULL, k-1, tt, c(0, Inf), data))
}

