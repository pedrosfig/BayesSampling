#' Bayes Linear Method for Categorical Data
#'
#' Creates the Bayes Linear Estimator for Categorical Data
#'
#' @param ys (k)-vector of sample proportions for each category.
#' @param n sample size.
#' @param N total size of the population.
#' @param m (k-1)-vector with the prior proportion of each strata. If NULL, sample proportion for each strata will be used (non-informative prior).
#' @param v vector with the prior variance - NÃO VAI TER
#' @param sigma vector with the prior estimate of variability (standard deviation) within each strata of the population. If NULL, sample variance of each strata will be used.
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
#' Estimator <- BLE_SSRS(ys,h,N,m,v,sigma)
#' Estimator
#' @export
BLE_SSRS <- function(ys, h, N, m=NULL,v=NULL,sigma=NULL){

  war_1 <- "parameter 'm' (prior proportions) not informed, sample proprotions used in estimations"
  war_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  war_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"


  # ys <- ys[-length(ys)] #retirando a observação da ultima categoria (se ys tiver tamanho k)
  # if(k == 1){stop("only 1 category defined, usar função da media amostral")}                  ###### ajustar

  k <- length(ys)


  # marker <- c(1)   #diz onde começam as obs de cada estrato
  # for(i in 2:H){
  #   marker <- c(marker, marker[i-1] + h[i-1])
  # }
  #
  # if (is.null(m)){
  #   warning(war_1)
  #   m <- c(mean(ys[1:h[1]]))
  #   for(i in 2:H-1){
  #     M <- mean(ys[marker[i]:(marker[i+1] - 1)])
  #     m <- c(m, M)
  #   }
  #   M <- mean(ys[marker[H] : length(ys)])
  #   m <- c(m, M)
  # }
  #
  #
  # if(is.null(sigma)){
  #   warning(war_2)
  #   s <- c(var(ys[1:h[1]]))
  #   for(i in 2:H-1){
  #     S <- var(ys[marker[i]:(marker[i+1] - 1)])
  #     s <- c(s, S)
  #   }
  #   S <- var(ys[marker[H] : length(ys)])
  #   s <- c(s, S)
  #   sigma <- sqrt(s)
  # }
  #
  # if(is.null(v)){
  #   warning(war_3)
  #   v <- c()
  #   for(i in 1:H){
  #     V <- 10^100 * m[i]
  #     v <- c(v, V)
  #   }
  # }
  #
  # for (i in 1:H) {
  #   if(v[i] < sigma[i]^2){
  #     stop("prior variance (parameter 'v') too small")
  #   }
  #
  # }


  a <- m
  v <- m*(1-m)
  xs <- diag(k)
  c <- v - sigma^2
  m_jj <- 1 -(sigma^2)/m
  Vs <-


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



  V_nots <- Vs*n/(N-n)

  R <- c*diag(H)

  return(BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots))



}
