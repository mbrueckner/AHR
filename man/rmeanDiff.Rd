% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmean.R
\name{rmeanDiff}
\alias{rmeanDiff}
\title{rmeanDiff}
\usage{
rmeanDiff(L, formula, data, rr.subset = rep(TRUE, nrow(data)))
}
\arguments{
\item{L}{time-limit specifying up to which time restricted mean will be calculated}

\item{formula}{an object of class '"formula"' specifying the conditional survival model}

\item{data}{data frame containing the variables in formula}

\item{rr.subset}{logical vector defining subset of observations to use for response rate estimation (default: use all observations)}
}
\value{
An object of class '"rmd"', i.e. a list containing:
 \item{L}{time limit, i.e. restricted mean up to time L is calculated}
 \item{rmean1}{restricted mean in group 1}
 \item{rmean2}{restricted mean in group 2}
 \item{rmean.diff}{estimated restricted mean difference}
 \item{var.rmean1}{an estimate of the asymptotic variance of the restricted mean in group 1}
 \item{var.rmean2}{an estimate of the asymptotic variance of the restricted mean in group 2}
 \item{var.rmean.diff}{an estimate of the asymptotic variance of the restricted mean difference}
 \item{Z.rmean}{the standardized test statistic for testing rmean.diff=0}
 \item{p.value}{p-value corresponding to Z.rmean}
}
\description{
Estimate difference of restricted mean survival based on (weighted) Kaplan-Meier estimates
of the survival functions in each group.
}
\examples{
T <- c(rexp(100, 1), rexp(100, 2))
C <- c(rexp(100, 1), rexp(100, 2))
Y <- pmin(T, C)
D <- T <= C
Z <- rep(c(0,1), c(100, 100))
fit <- rmeanDiff(2, formula=Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z))
}
