#' creates vector of 1's to be used in the estimators
#' @param y sample matrix
#' @return vector of 1's with size equal to the number of observations in the sample
#'
create1 <- function(y){
  if(is.vector(y)){
    vect1 <- rep(1,length(y))
  }
  else if(is.list(y)){
    vect1 <- rep(1,dim(y)[1])
  }
  else if(is.matrix(y)){
    vect1 <- rep(1,length(y[,1]))
  }
  else {
    stop("incorrect number of dimensions")
  }
  return(vect1)
  }
