#' Bayes Linear Method for Categorical Data
#'
#' Creates the Bayes Linear Estimator for Categorical Data
#'
#' @param ys k-vector of sample proportion for each category.
#' @param n sample size.
#' @param N total size of the population.
#' @param m k-vector with the prior proportion of each strata. If \code{NULL}, sample proportion for each strata will be used (non-informative prior).
#' @param rho matrix with the prior correlation coefficients between two different units within categories. It must be a symmetric square matrix of dimension k (or k-1). If \code{NULL}, non-informative prior will be used.
#'
#' @return A list containing the following components: \itemize{
#' \item \code{est.prop} - BLE for the sample proportion of each category
#' \item \code{Vest.prop} - Variance associated with the above
#' \item \code{Vs.Matrix} - Vs matrix, as defined by the BLE method (should be a positive-definite matrix)
#' \item \code{R.Matrix} - R matrix, as defined by the BLE method (should be a positive-definite matrix)
#' }
#'
#' @source \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400111886}
#' @references Gonçalves, K.C.M, Moura, F.A.S and  Migon, H.S.(2014). Bayes Linear Estimation for Finite Population with emphasis on categorical data. Survey Methodology, 40, 15-28.
#'
#' @examples
#' # 2 categories
#' ys <- c(0.2614, 0.7386)
#' n <- 153
#' N <- 15288
#' m <- c(0.7, 0.3)
#' rho <- matrix(0.1, 1)
#'
#' Estimator <- BLE_Categorical(ys,n,N,m,rho)
#' Estimator
#' 
#' ys <- c(0.2614, 0.7386)
#' n <- 153
#' N <- 15288
#' m <- c(0.7, 0.3)
#' rho <- matrix(0.5, 1)
#'
#' Estimator <- BLE_Categorical(ys,n,N,m,rho)
#' Estimator
#' 
#' 
#' # 3 categories
#' ys <- c(0.2, 0.5, 0.3)
#' n <- 100
#' N <- 10000
#' m <- c(0.4, 0.1, 0.5)
#' mat <- c(0.4, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.6)
#' rho <- matrix(mat, 3, 3)
#'
#' @export
BLE_Categorical <- function(ys, n, N, m=NULL, rho=NULL){

  mes_1 <- "parameter 'm' (prior proportions) not informed, sample proportions used in estimations"
  mes_2 <- "parameter 'rho' not informed, non informative prior correlation coefficients used in estimations"
  
  k <- length(ys)
  if( k == 1 ){stop("only 1 category defined")}
  
  if( prod(ys >= 0) != 1 ){stop("all sample proportions must be non-negative numbers")}
  if( sum(ys) != 1 ){stop("sum of sample proportions should be 1")}
  
  if( is.null(m) ){
    message(mes_1)
    m <- ys
  }
  
  if( length(m) != length(ys) ){stop("length of parameters 'ys' and 'm' must coincide")}
  if( prod(m > 0) != 1 ){stop("all prior proportions must be positive numbers")}
  if( sum(m) != 1 ){stop("sum of prior proportions should be 1")}
  
  rho_informed <- 1
  if( is.null(rho) ){
    rho_informed <- 0
    message(mes_2)
    rho <- diag(x = 1-1e-10, nrow=k)
    }
  
  if( ! is.symmetric.matrix(rho) | dim(rho)[1] < k-1 ){stop("'rho' must be a symmetric square matrix of dimension k")}
  if( max(abs(rho)) >= 1 ){stop("all values in 'rho' must be between -1 and 1")}
  
  ys <- ys[-k]
  m <- m[-k]
  
  
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
  
  if( rho_informed == 1 ){
    if( ! is.symmetric.matrix(R) ){stop("'R' must be a symmetric matrix. Review parameter 'rho'")}
    if( ! is.positive.definite(R, tol=1e-15) ){warning("'R' should be a positive-definite matrix. Possible problem with parameter 'rho'")}
  }
  
  
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
  
  if( rho_informed == 1 ){
    if( ! is.symmetric.matrix(Vs) ){stop("'Vs' must be a symmetric matrix. Review parameter 'rho'")}
    if( ! is.positive.definite(Vs, tol=1e-15) ){warning("'Vs' should be a positive-definite matrix. Possible problem with parameter 'rho'")}
  }
  

  V_nots <- Vs*n/(N-n)
  
  C_inv <- ginv(R) + ginv(Vs)
  C <- ginv(C_inv)
  Beta <- C%*%(ginv(Vs)%*%ys + ginv(R)%*%a)
  
  
  p_aux <- (n*ys + (N-n)*Beta)/N
  p_k <- 1 - sum(p_aux)
  p <- c(p_aux,p_k)
  
  V_aux <- (V_nots + C) * ((N-n)/N)^2
  V_k <- sum(V_aux)
  Cov_k <- c
  for (i in 1:k-1) {
    Cov_k[i] <- -sum(V_aux[i,])
  }
  V_p <- rbind(cbind(V_aux, Cov_k[]), c(Cov_k, V_k))
  
  if( prod(diag(V_p) > 0) != 1 ){warning("'Vest.prop' should have only positive diagonal values. Review prior parameters.")}
  
  return(list(est.prop = p, Vest.prop = V_p, Vs.Matrix = Vs, R.Matrix = R))
  
}

