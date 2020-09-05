#' Simple Random Sample BLE
#'
#' Creates the Bayes Linear Estimator for the Simple Random Sampling design (without replacement)
#'
#' @param ys vector of sample observations or sample mean ('sigma' and 'n' parameters will be required in this case).
#' @param N total size of the population.
#' @param n sample size. Necessary only if 'ys' represents sample mean (will not be used otherwise).
#' @param m prior mean. If NULL, sample mean will be used (non-informative prior).
#' @param v prior variance of an element from the population (bigger than sigma^2). If NULL, it will tend to infinity (non-informative prior).
#' @param sigma prior estimate of variability (standard deviation) within the population. If NULL, sample variance will be used.
#'
#' @return A list containing the following components: \itemize{
#' \item \code{est.beta} - BLE of Beta (BLE for every individual)
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
#' ys <- c(5,6,8)
#' m <- 6
#' v <- 5
#' sigma <- 1
#' N <- 5
#'
#' Estimator <- BLE_SRS(ys,N,m,v,sigma)
#' Estimator
#' @export
BLE_SRS <- function(ys, N, n=NULL, m=NULL, v=NULL, sigma=NULL){

  war_1 <- "parameter 'm' (prior mean) not informed, sample mean used in estimations"
  war_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  war_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"

  if (is.null(m)){
    warning(war_1)
    m <- mean(ys)
  }


  if(length(ys)==1){
    if( (is.null(sigma)) | is.null(n) ){
      stop("ys of length 1 requires not null parameters 'sigma' and 'n'")
    }
    ys <- rep(ys, n)
  }


  if(is.null(sigma)){
    warning(war_2)
    sigma <- sqrt(var(ys))
  }

  if(is.null(v)){
    warning(war_3)
    v <- 10^100 * mean(ys)}

  if(v < sigma^2){
    stop("prior variance (parameter 'v') too small")
  }


  xs <- create1(ys)
  a <- m
  Vs <- diag(xs)*(sigma^2)
  c <- v - sigma^2
  R <- c

  x_nots <- rep(1,N - length(ys))
  V_nots <- diag(x_nots)*(sigma^2)

  return(BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots))
}




#' Stratified Simple Random Sample BLE
#'
#' Creates the Bayes Linear Estimator for the Stratified Simple Random Sampling design (without replacement)
#' @param ys vector of sample observations.
#' @param h vector with number of observations in each strata.
#' @param N vector with the total size of each strata.
#' @param m vector with the prior mean of each strata. If NULL, sample mean for each strata will be used (non-informative prior).
#' @param v vector with the prior variance of an element from each strata (bigger than sigma^2 for each strata). If NULL, it will tend to infinity (non-informative prior).
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

  war_1 <- "parameter 'm' (prior mean) not informed, sample mean used in estimations"
  war_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  war_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"



  H <- length(h)
  if(H == 1){stop("only 1 strata defined, try using the BLE_SRS() function")}


  marker <- c(1)   #diz onde começam as obs de cada estrato
  for(i in 2:H){
    marker <- c(marker, marker[i-1] + h[i-1])
  }

  if (is.null(m)){
    warning(war_1)
    m <- c(mean(ys[1:h[1]]))
    for(i in 2:H-1){
      M <- mean(ys[marker[i]:(marker[i+1] - 1)])
      m <- c(m, M)
    }
    M <- mean(ys[marker[H] : length(ys)])
    m <- c(m, M)
  }


  if(is.null(sigma)){
    warning(war_2)
    s <- c(var(ys[1:h[1]]))
    for(i in 2:H-1){
      S <- var(ys[marker[i]:(marker[i+1] - 1)])
      s <- c(s, S)
    }
    S <- var(ys[marker[H] : length(ys)])
    s <- c(s, S)
    sigma <- sqrt(s)
  }

  if(is.null(v)){
    warning(war_3)
    v <- c()
    for(i in 1:H){
      V <- 10^100 * m[i]
      v <- c(v, V)
    }
  }

  for (i in 1:H) {
    if(v[i] < sigma[i]^2){
      stop("prior variance (parameter 'v') too small")
    }

  }


  aux <- rep(1, h[1])
  for(i in 2:H){
    zero <- rep(0, sum(h))
    one <- rep(1, h[i])
    aux <- c(aux, zero, one)
  }
  xs <- matrix(aux,nrow = sum(h),ncol=H)


  out <- N-h
  aux_out <- rep(1, out[1])
  for(i in 2:H){
    zero <- rep(0, sum(out))
    one <- rep(1, out[i])
    aux_out <- c(aux_out, zero, one)
  }
  x_nots <- matrix(aux_out,nrow = sum(out),ncol=H)


  Vs <- diag(h[1])*(sigma[1])^2
  for (i in 2:H) {
    V <- diag(h[i])*(sigma[i])^2
    Vs <- bdiag(Vs,V)
  }
  Vs <- as.matrix(Vs)


  k <- N[1] - h[1]
  V_nots <- diag(k)*(sigma[1])^2
  for (i in 2:H) {
    k <- N[i] - h[i]
    V <- diag(k)*(sigma[i])^2
    V_nots <- bdiag(V_nots,V)
  }
  V_nots <- as.matrix(V_nots)


  a <- m
  c <- v - sigma^2
  R <- c*diag(H)

  return(BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots))



}




#' Ratio BLE
#'
#' Creates the Bayes Linear Estimator for the Ratio "estimator"
#'
#' @param ys vector of sample observations.
#' @param xs vector with values for the auxiliary variable of the elements in the sample.
#' @param x_nots vector with values for the auxiliary variable of the elements not in the sample.
#' @param m prior mean for the ratio between Y and X. If NULL, \code{mean(ys)/mean(xs)} will be used (non-informative prior).
#' @param v prior variance of the ratio between Y and X (bigger than sigma^2). If NULL, it will tend to infinity (non-informative prior).
#' @param sigma prior estimate of variability (standard deviation) of the ratio within the population. If NULL, sample variance of the ratio will be used.
#'
#' @return A list containing the following components: \itemize{
#' \item \code{est.beta} - BLE of Beta
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
#' ys <- c(10,8,6)
#' xs <- c(5,4,3)
#' x_nots <- c(1,20)
#' m <- 2
#' v <- 10
#' sigma <- 2
#'
#' Estimator <- BLE_Ratio(ys,xs,x_nots,m,v,sigma)
#' Estimator
#' @export
BLE_Ratio <- function(ys,xs, x_nots,m=NULL,v=NULL,sigma=NULL){

  war_1 <- "parameter 'm' (prior mean) not informed, sample mean used in estimations"
  war_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  war_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"


  z <- ys/xs

  if (is.null(m)){
    warning(war_1)
    m <- mean(ys)/mean(xs)
  }

  if(is.null(sigma)){
    warning(war_2)
    sigma <- sqrt(var(z))
  }

  if(is.null(v)){
    warning(war_3)
    v <- 10^100 * mean(ys)}

  if(v < sigma^2){
    stop("prior variance (parameter 'v') too small")
  }


  Vs <- diag(xs)*(sigma^2)
  V_nots <- diag(x_nots)*(sigma^2)

  a <- m
  c <- v - sigma^2
  R <- c

  return(BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots))



}
















