#' UpdateConst
#'
#' updates the vector const after inputting sample y
#'
#'
#' @name UpdateConst
#' @param y A matrix, the new sample to be input.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @return The updated vector which specifies all the constant parameters.





UpdateConst <- function (y,const)
{
  const[1] = const[1] + dim(y)[1]
  const
}

