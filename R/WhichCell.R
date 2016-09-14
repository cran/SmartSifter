#' WhichCell
#'
#' returns the index of the cell to which the new sample belongs
#'
#'
#' @name WhichCell
#' @param y A one-row matrix, the new sample to be input.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @return The row index of the discrete class to which the new sample belongs.


WhichCell <- function (y, param, smart, const)
{
  n = const[2]
  tmpsmart = as.matrix(smart[,1:n])
  tmpy = y[,1:n]

  numOfCell=dim(smart)[1]

  index = numOfCell + 1

  for (i in 1:numOfCell){
    if (tmpsmart[i,]==tmpy) {index=i}
  }

  index
}
