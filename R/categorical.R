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
BLE_Categorical <- function(ys, h, N, m=NULL, rho=NULL){

  war_1 <- "parameter 'm' (prior proportions) not informed, sample proprotions used in estimations"
  war_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  

  k <- length(ys)
  if(k == 1){stop("only 1 category defined")}
  
  if(sum(ys) != 1){stop("sum of sample proportions should be 1")}
  
  ys <- ys[-k]      #retirando a observação da ultima categoria
  m <- m[-k]        #retirando a priori da ultima categoria
  
  
  
  if( ! is.symmetric.matrix(rho) ){stop("rho must be a symmetric square matrix")}
  
  
  m_ij <- matrix(0, nrow = k, ncol = k)
  
  for (i in 1:(k-1)) {
    for (j in 1:(k-1)) {
      m_ij[i,j] <- (m[i]*m[j] + rho[j,i]*sqrt( m[i]*(1-m[i])*m[j]*(1-m[j]) ))/m[j]
    }
  }
  
  
  

  a <- m
  v <- m*(1-m)    # ver se é necessário !!!!
  c <- m * (diag(m_ij) - m)
#  sigma <- sqrt( m*(1-diag(m_ij)) )
  sigma <- sqrt(v - c)

  xs <- diag(k-1)
  
  
  # Construção do Vs
  
  Vs_d <- diag(sigma^2, nrow = k-1)           # elements in the diagonal
  
  Vs_out <- matrix(0, nrow = k-1, ncol = k-1)
  for (i in 1:(k-1)) {                        # elements outside the diagonal
    for (j in 1:(k-1)) {
      Vs_out[i,j] <- (-m[i])*m_ij[j,i]
    }
  }
  Vs_out <- Vs_out * (1 - diag(nrow = k-1))
  
  Vs <- (1/n)*(Vs_d + Vs_out)
  #  if( ! is.symmetric.matrix(Vs) ){stop("Vs must be a symmetric matrix. Review parameter 'rho'")}   - proximo teste ja vai dar erro
  if( ! is.positive.definite(Vs) ){stop("Vs must be a positive-definite matrix. Review parameter 'rho'")}
  

  V_nots <- Vs*n/(N-n)  # checar
  
  
  
  # Construção do R 
  
  R_d <- diag(c, nrow = k-1)                  # elements in the diagonal
  
  R_out <- matrix(0, nrow = k-1, ncol = k-1)
  for (i in 1:(k-1)) {                        # elements outside the diagonal
    for (j in 1:(k-1)) {
      R_out[i,j] <- m[i] * (m_ij[j,i] - m[j])
    }
  }
  R_out <- R_out * (1 - diag(nrow = k-1))
  
  R <- (1/n)*(R_d + R_out)
  #  if( ! is.symmetric.matrix(R) ){stop("R must be a symmetric matrix. Review parameter 'rho'")}   - proximo teste ja vai dar erro
  if( ! is.positive.definite(R) ){stop("R must be a positive-definite matrix. Review parameter 'rho'")}
  
  
  

#
#   out <- N-h
#   aux_out <- rep(1, out[1])
#   for(i in 2:H){
#     zero <- rep(0, sum(out))
#     one <- rep(1, out[i])
#     aux_out <- c(aux_out, zero, one)
#   }
#   x_nots <- matrix(aux_out,nrow = sum(out),ncol=H)


  # Vs <- diag(h[1])*(sigma[1])^2
  # for (i in 2:H) {
  #   V <- diag(h[i])*(sigma[i])^2
  #   Vs <- bdiag(Vs,V)
  # }
  # Vs <- as.matrix(Vs)





  return(BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots))



}
