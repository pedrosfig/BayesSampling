% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categorical.R
\name{BLE_Categorical}
\alias{BLE_Categorical}
\title{Bayes Linear Method for Categorical Data}
\source{
\url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400111886}
}
\usage{
BLE_Categorical(ys, n, N, m = NULL, rho = NULL)
}
\arguments{
\item{ys}{k-vector of sample proportion for each category.}

\item{n}{sample size.}

\item{N}{total size of the population.}

\item{m}{k-vector with the prior proportion of each strata. If \code{NULL}, sample proportion for each strata will be used (non-informative prior).}

\item{rho}{matrix with the prior correlation coefficients between two different units within categories. It must be a symmetric square matrix of dimension k (or k-1). If \code{NULL}, non-informative prior will be used.}
}
\value{
A list containing the following components: \itemize{
\item \code{est.prop} - BLE for the sample proportion of each category
\item \code{Vest.prop} - Variance associated with the above
\item \code{Vs.Matrix} - Vs matrix, as defined by the BLE method (should be a positive-definite matrix)
\item \code{R.Matrix} - R matrix, as defined by the BLE method (should be a positive-definite matrix)
}
}
\description{
Creates the Bayes Linear Estimator for Categorical Data
}
\examples{
# 2 categories
ys <- c(0.2614, 0.7386)
n <- 153
N <- 15288
m <- c(0.7, 0.3)
rho <- matrix(0.1, 1)

Estimator <- BLE_Categorical(ys,n,N,m,rho)
Estimator

ys <- c(0.2614, 0.7386)
n <- 153
N <- 15288
m <- c(0.7, 0.3)
rho <- matrix(0.5, 1)

Estimator <- BLE_Categorical(ys,n,N,m,rho)
Estimator


# 3 categories
ys <- c(0.2, 0.5, 0.3)
n <- 100
N <- 10000
m <- c(0.4, 0.1, 0.5)
mat <- c(0.4, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.6)
rho <- matrix(mat, 3, 3)

}
\references{
Gonçalves, K.C.M, Moura, F.A.S and  Migon, H.S.(2014). Bayes Linear Estimation for Finite Population with emphasis on categorical data. Survey Methodology, 40, 15-28.
}
