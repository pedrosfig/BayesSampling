# Cap 3.2


#' calculates the C factor
#'
#'@inheritParams BLE_Reg

C <- function(ys,xs,R,Vs) {
  c1 <- ginv(R)
  c2 <- t(xs)%*%ginv(Vs)%*%xs
  C_minus <- c1 + c2
  C_result <- ginv(C_minus)
  return(C_result)
}


#' calculates the BLE for Beta
#'
#'@inheritParams BLE_Reg

E_beta <- function(ys,xs,a,R,Vs) {
  c_beta <- C(ys,xs,R,Vs)
  p1 <- t(xs)%*%ginv(Vs)%*%ys
  p2 <- ginv(R)%*%a
  res <- c_beta%*%(p1 + p2)
  return(res)
}


#' calculates the risk matrix associated with the BLE for Beta
#'
#'@inheritParams BLE_Reg

V_beta <- function(ys,xs,R,Vs){
  v_beta <- C(ys,xs,R,Vs)
  return(v_beta)
}


#'calculates the BLE for the individuals not in the sample
#'@inheritParams BLE_Reg

E_theta_Reg <- function(ys,xs,a,R,Vs,x_nots) {
  c_theta <- C(ys,xs,R,Vs)
  p1 <- t(xs)%*%ginv(Vs)%*%ys
  p2 <- ginv(R)%*%a
  res <- x_nots%*%c_theta%*%(p1 + p2)
  return(res)
}


#'calculates the risk matrix associated with the BLE for the individuals not in the sample
#'@inheritParams BLE_Reg

V_theta_Reg <- function(ys,xs,R,Vs,x_nots,V_nots) {
  c_theta <- C(ys,xs,R,Vs)
  p1 <- x_nots%*%c_theta%*%t(x_nots)
  p2 <- V_nots
  res <- p1 + p2
  return(res)
}


#'calculates BLE for the total T
#'@inheritParams BLE_Reg
#'

T_Reg <- function(ys,xs,a,R,Vs,x_nots) {
  one_s <- create1(ys)
  parc1 <- t(one_s)%*%ys
  esp <- E_theta_Reg(ys,xs,a,R,Vs,x_nots)
  one_nots <- create1(esp)
  parc2 <- t(one_nots)%*%esp
  return(parc1 + parc2)
}


#'calculates risk matrix associated with the BLE for for the total T
#'@inheritParams BLE_Reg
#'

VT_Reg <- function(ys,xs,a,R,Vs,x_nots,V_nots) {
  v_theta <- V_theta_Reg(ys,xs,R,Vs,x_nots,V_nots)
  esp <- E_theta_Reg(ys,xs,a,R,Vs,x_nots)
  one_nots <- create1(esp)
  res <- t(one_nots)%*%v_theta%*%one_nots
  return(res)
}


#' General BLE case
#'
#' Calculates the Bayes Linear Estimator for Regression models (general case)
#'@param ys response variable of the sample
#'@param xs explicative variable of the sample
#'@param a vector of means from Beta
#'@param R covariance matrix of Beta
#'@param Vs covariance of sample errors
#'@param x_nots values of X for the individuals not in the sample
#'@param V_nots covariance matrix of the individuals not in the sample
#'
#' @return A list containing the following components: \itemize{
#' \item est.beta - BLE of Beta
#' \item Vest.beta - Variance associated with the above
#' \item est.mean - BLE of each individual not in the sample
#' \item Vest.mean - Covariance matrix associated with the above
#' \item est.tot - BLE for the total
#' \item Vest.tot - Variance associated with the above
#' }
#'
#' @source \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400111886}
#' @references Gonçalves, K.C.M, Moura, F.A.S and  Migon, H.S.(2014). Bayes Linear Estimation for Finite Population with emphasis on categorical data. Survey Methodology, 40, 15-28.
#'
#' @examples
#' xs <- matrix(c(1,1,1,1,2,3,5,0),nrow=4,ncol=2)
#' ys <- c(12,17,28,2)
#' x_nots <- matrix(c(1,1,1,0,1,4),nrow=3,ncol=2)
#' a <- c(1.5,6)
#' R <- matrix(c(10,2,2,10),nrow=2,ncol=2)
#' Vs <- diag(c(1,1,1,1))
#' V_nots <- diag(c(1,1,1))
#'
#' Estimator <- BLE_Reg(ys,xs,a,R,Vs,x_nots,V_nots)
#' Estimator
#'
#' @export
BLE_Reg <- function(ys,xs,a,R,Vs,x_nots,V_nots){
  beta <- as.data.frame(E_beta(ys,xs,a,R,Vs))
  colnames(beta) = c("Beta")
  var_beta <- as.data.frame(V_beta(ys,xs,R,Vs))

  y_nots <- as.data.frame(E_theta_Reg(ys,xs,a,R,Vs,x_nots))
  colnames(y_nots) = c("y_nots")
  var_y_nots <- as.data.frame(V_theta_Reg(ys,xs,R,Vs,x_nots,V_nots))

  total <- T_Reg(ys,xs,a,R,Vs,x_nots)[,]
  var_total <- VT_Reg(ys,xs,a,R,Vs,x_nots,V_nots)[,]

  return(list(est.beta = beta, Vest.beta = var_beta,est.mean = y_nots, Vest.mean = var_y_nots, est.tot = total, Vest.tot = var_total))
  }

