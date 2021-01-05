#' Simple Random Sample BLE
#'
#' Creates the Bayes Linear Estimator for the Simple Random Sampling design (without replacement)
#'
#' @param ys vector of sample observations or sample mean (\code{sigma} and \code{n} parameters will be required in this case).
#' @param N total size of the population.
#' @param m prior mean. If \code{NULL}, sample mean will be used (non-informative prior).
#' @param v prior variance of an element from the population (bigger than \code{sigma^2}). If \code{NULL}, it will tend to infinity (non-informative prior).
#' @param sigma prior estimate of variability (standard deviation) within the population. If \code{NULL}, sample variance will be used.
#' @param n sample size. Necessary only if \code{ys} represent sample mean (will not be used otherwise).
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
#' N <- 5
#' m <- 6
#' v <- 5
#' sigma <- 1
#'
#' Estimator <- BLE_SRS(ys, N, m, v, sigma)
#' Estimator
#'
#'
#' # Same example but informing sample mean and sample size instead of sample observations
#' ys <- mean(c(5,6,8))
#' N <- 5
#' n <- 3
#' m <- 6
#' v <- 5
#' sigma <- 1
#'
#' Estimator <- BLE_SRS(ys, N, m, v, sigma, n)
#' Estimator
#'
#' @export
BLE_SRS <- function(ys, N, m=NULL, v=NULL, sigma=NULL, n=NULL){

  mes_1 <- "parameter 'm' (prior mean) not informed, sample mean used in estimations"
  mes_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  mes_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"
  mes_4 <- "sample mean informed instead of sample observations, parameters 'n' and 'sigma' will be necessary"


  if(length(ys)==1){
    message(mes_4)
    if( (is.null(sigma)) | is.null(n) ){
      stop("ys of length 1 requires not null parameters 'sigma' and 'n'")
    }
    ys <- rep(ys, n)
  }


  if (is.null(m)){
    message(mes_1)
    m <- mean(ys)
  }

  if(is.null(sigma)){
    message(mes_2)
    sigma <- sqrt(var(ys))
  }

  if(is.null(v)){
    message(mes_3)
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
#' @param ys vector of sample observations or sample mean for each strata (\code{sigma} parameter will be required in this case).
#' @param h vector with number of observations in each strata.
#' @param N vector with the total size of each strata.
#' @param m vector with the prior mean of each strata. If \code{NULL}, sample mean for each strata will be used (non-informative prior).
#' @param v vector with the prior variance of an element from each strata (bigger than \code{sigma^2} for each strata). If \code{NULL}, it will tend to infinity (non-informative prior).
#' @param sigma vector with the prior estimate of variability (standard deviation) within each strata of the population. If \code{NULL}, sample variance of each strata will be used.
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
#' N <- c(5,5,3)
#' m <- c(0,9,8)
#' v <- c(3,8,1)
#' sigma <- c(1,2,0.5)
#'
#' Estimator <- BLE_SSRS(ys, h, N, m, v, sigma)
#' Estimator
#'
#'
#' # Same example but informing sample means instead of sample observations
#' y1 <- mean(c(2,-1,1.5))
#' y2 <- mean(c(6,10))
#' y3 <- mean(c(8,8))
#' ys <- c(y1, y2, y3)
#' h <- c(3,2,2)
#' N <- c(5,5,3)
#' m <- c(0,9,8)
#' v <- c(3,8,1)
#' sigma <- c(1,2,0.5)
#'
#' Estimator <- BLE_SSRS(ys, h, N, m, v, sigma)
#' Estimator
#'
#' @export
BLE_SSRS <- function(ys, h, N, m=NULL, v=NULL, sigma=NULL){

  mes_1 <- "parameter 'm' (prior mean) not informed, sample mean used in estimations"
  mes_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  mes_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"
  mes_4 <- "sample means informed instead of sample observations, parameter 'sigma' will be necessary"


  H <- length(h)
  if(H == 1){stop("only 1 strata defined, try using the BLE_SRS() function")}


  if(length(ys)!=sum(h)){
    if(length(ys)!=length(h)){
      stop("length of 'ys' incompatable with 'h'")
    }
    else{                # length(ys)==length(h)
      message(mes_4)
      if(is.null(sigma)){
        stop("not null parameter 'sigma' required")
      }
      ys <- rep(ys, h)
    }
  }


  marker <- c(1)   # marks where the observations of each strata begin
  for(i in 2:H){
    marker <- c(marker, marker[i-1] + h[i-1])
  }


  if (is.null(m)){
    message(mes_1)
    m <- c(mean(ys[1:h[1]]))
    for(i in 2:H-1){
      M <- mean(ys[marker[i]:(marker[i+1] - 1)])
      m <- c(m, M)
    }
    M <- mean(ys[marker[H] : length(ys)])
    m <- c(m, M)
  }

  if(is.null(sigma)){
    message(mes_2)
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
    message(mes_3)
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
#' @param ys vector of sample observations or sample mean (\code{sigma} and \code{n} parameters will be required in this case).
#' @param xs vector with values for the auxiliary variable of the elements in the sample or sample mean.
#' @param x_nots vector with values for the auxiliary variable of the elements not in the sample.
#' @param m prior mean for the ratio between Y and X. If \code{NULL}, \code{mean(ys)/mean(xs)} will be used (non-informative prior).
#' @param v prior variance of the ratio between Y and X (bigger than \code{sigma^2}). If \code{NULL}, it will tend to infinity (non-informative prior).
#' @param sigma prior estimate of variability (standard deviation) of the ratio within the population. If \code{NULL}, sample variance of the ratio will be used.
#' @param n sample size. Necessary only if \code{ys} and \code{xs} represent sample means (will not be used otherwise).
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
#' xs <- c(5,4,3.1)
#' x_nots <- c(1,20,13,15,-5)
#' m <- 2.5
#' v <- 10
#' sigma <- 2
#'
#' Estimator <- BLE_Ratio(ys, xs, x_nots, m, v, sigma)
#' Estimator
#'
#'
#' # Same example but informing sample means and sample size instead of sample observations
#' ys <- mean(c(10,8,6))
#' xs <- mean(c(5,4,3.1))
#' n <- 3
#' x_nots <- c(1,20,13,15,-5)
#' m <- 2.5
#' v <- 10
#' sigma <- 2
#'
#' Estimator <- BLE_Ratio(ys, xs, x_nots, m, v, sigma, n)
#' Estimator
#'
#' @export
BLE_Ratio <- function(ys, xs, x_nots, m=NULL, v=NULL, sigma=NULL, n=NULL){

  mes_1 <- "parameter 'm' (prior mean) not informed, sample mean used in estimations"
  mes_2 <- "parameter 'sigma' (prior variability) not informed, sample variance used in estimations"
  mes_3 <- "parameter 'v' (prior variance of an element) not informed, (10^100 * mean(ys)) used in estimations (non-informative prior)"
  mes_4 <- "sample means informed instead of sample observations, parameters 'n' and 'sigma' will be necessary"


  if(length(ys) != length(xs)){
    stop("dimensions of ys and xs are different")
  }


  if(length(ys)==1){
    message(mes_4)
    if( (is.null(sigma)) | is.null(n) ){
      stop("ys of length 1 requires not null parameters 'sigma' and 'n'")
    }
    ys <- rep(ys, n)
    xs <- rep(xs, n)
  }


  z <- ys/xs


  if (is.null(m)){
    message(mes_1)
    m <- mean(ys)/mean(xs)
  }

  if(is.null(sigma)){
    message(mes_2)
    sigma <- sqrt(var(z))
  }

  if(is.null(v)){
    message(mes_3)
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
















