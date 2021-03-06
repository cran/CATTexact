% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CATTexact.R
\name{catt_exact}
\alias{catt_exact}
\title{Conditional exact Cochran-Armitage trend test}
\usage{
catt_exact(dose.ratings, totals, cases)
}
\arguments{
\item{dose.ratings}{A vector of dose ratings, the i-th entry corresponds to the dose-rating of the i-th group. This vector must be strictly monotonically increasing}

\item{totals}{The vector of total individuals per group, the i-th entry corresponds to the total number of individuals in the i-th group.}

\item{cases}{The vector of incidences per groups, the i-th entry corresponds to the number of incidences in the i-th group.}
}
\value{
A list containing the value of the Cochran-Armitage Trend Test Statistic, its exact and asymptotic p-value.
}
\description{
\code{catt_exact} calculates the Cochran-Armitage trend test statistic (Cochran (1954), Armitage (1955)) and the one-sided p-value for the corresponding conditional exact test.
The conditional exact test has been established by Williams (1988). The computation of its p-value is performed using an algorithm following an idea by Mehta, et al. (1992).
}
\examples{
d <- c(1,2,3,4)
n <- rep(20,4)
r <- c(1,4,3,8)

catt_exact(d, n, r)

}
\references{
Armitage, P. Tests for linear trends in proportions and frequencies. \emph{Biometrics}, 11 (1955): 375-386.

Cochran, W. G. Some methods for strengthening the common \eqn{\chi^2} tests, \emph{Biometrics}. 10 (1954): 417-451.

Mehta, C. R., Nitin P., and Pralay S. Exact stratified linear rank tests for ordered categorical and binary data. \emph{Journal of Computational and Graphical Statistics}, 1 (1992): 21-40.

Portier, C., and Hoel D. Type 1 error of trend tests in proportions and the design of cancer screens. \emph{Communications in Statistics-Theory and Methods}, 13 (1984): 1-14.

Williams, D. A. Tests for differences between several small proportions. \emph{Applied Statistics}, 37 (1988): 421-434.
}
