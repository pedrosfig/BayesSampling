#' Bayes Linear Method for Categorical Data
#'
#' Creates the Bayes Linear Estimator for Categorical Data
#'
#' @param ys k-vector of sample proportion for each category.
#' @param n sample size.
#' @param N total size of the population.
#' @param m k-vector with the prior proportion of each strata. If NULL, sample proportion for each strata will be used (non-informative prior).
#' @param rho matrix with the prior correlation coefficients between two different units within categories. It must be a symmetric square matrix of dimension k.
#'
#' @return A list containing the following components: \itemize{
#' \item \code{est.beta} - BLE of Beta (BLE for the individuals in each strata)
#' \item \code{Vest.beta} - Variance associated with the above
#' \item \code{est.mean} - BLE for each individual not in the sample
#' \item \code{Vest.mean} - Covariance matrix associated with the above
#' \item \code{est.tot} - BLE for the total
#' \item \code{Vest.tot} - Variance associated with the above
#' }
#'
#' @source \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400111886}
#' @references Gonçalves, K.C.M, Moura, F.A.S and  Migon, H.S.(2014). Bayes Linear Estimation for Finite Population with emphasis on categorical data. Survey Methodology, 40, 15-28.
#'
#' @examples
#' ys <- c(2,-1,1.5, 6,10, 8,8)
#' h <- c(3,2,2)
#' m <- c(0,9,8)
#' v <- c(3,8,1)
#' sigma <- c(1,2,0.5)
#' N <- c(5,5,3)
#'
#' Estimator <- BLE_SSRS(ys,h,N,m,rho)
#' Estimator
#' @export
BLE_Categorical <- function(ys, n, N, m=NULL, rho=NULL){

  mes_1 <- "parameter 'm' (prior proportions) not informed, sample proprotions used in estimations"
   

  k <- length(ys)
  if(k == 1){stop("only 1 category defined")}
  
  if(prod(ys >= 0) != 1){stop("all sample proportions must be non-negative numbers")}
  if(sum(ys) != 1){stop("sum of sample proportions should be 1")}
  
  ys <- ys[-k]
  m <- m[-k]
  
  
  if( ! is.symmetric.matrix(rho) ){stop("rho must be a symmetric square matrix")}
  
  
  # Matrix 'm_ij'
  
  m_ij <- matrix(0, nrow = k-1, ncol = k-1)
  
  for (i in 1:(k-1)) {
    for (j in 1:(k-1)) {
      m_ij[i,j] <- (m[i]*m[j] + rho[j,i]*sqrt( m[i]*(1-m[i])*m[j]*(1-m[j]) ))/m[j]
    }
  }
  
  

  a <- m
  v <- m*(1-m)
  c <- m * (diag(m_ij) - m)
  sigma <- sqrt(v - c)

  xs <- diag(k-1)
  x_nots <- diag(k-1)
  
  
  
  # Matrix 'R' 
  
  R_d <- diag(c, nrow = k-1)                    # elements in the diagonal
  
  R_out <- matrix(0, nrow = k-1, ncol = k-1)
  for (i in 1:(k-1)) {                        
    for (j in 1:(k-1)) {
      R_out[i,j] <- m[i] * (m_ij[j,i] - m[j])   # elements outside the diagonal
    }
  }
  R_out <- R_out * (1 - diag(nrow = k-1))
  
  R <- R_d + R_out
  
  if( ! is.symmetric.matrix(R) ){stop("R must be a symmetric matrix. Review parameter 'rho'")}
  if( ! is.positive.definite(R) ){stop("R must be a positive-definite matrix. Review parameter 'rho'")}
  
  
  
  # Matrix 'Vs'
  
  Vs_d <- diag(sigma^2, nrow = k-1)             # elements in the diagonal
  
  Vs_out <- matrix(0, nrow = k-1, ncol = k-1)
  for (i in 1:(k-1)) {                        
    for (j in 1:(k-1)) {
      Vs_out[i,j] <- (-m[i])*m_ij[j,i]          # elements outside the diagonal
    }
  }
  Vs_out <- Vs_out * (1 - diag(nrow = k-1))
  
  Vs <- (1/n)*(Vs_d + Vs_out)
  
  if( ! is.symmetric.matrix(Vs) ){stop("Vs must be a symmetric matrix. Review parameter 'rho'")}
  if( ! is.positive.definite(Vs) ){stop("Vs must be a positive-definite matrix. Review parameter 'rho'")}
  

  V_nots <- Vs*n/(N-n)
  
  
  
  return(BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots))



}
