#' InputSample
#'
#' updates parameters after inputting sample (can be more than one)
#'
#'
#' @name InputSample
#' @param y A matrix, the new sample to be input.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @return The matrix which stores updated parameters over the continuous domain.



InputSample <- function (y, param = TRUE, smart, const, initial)
{
  numOfSample = dim(y)[1]

  for (i in 1:numOfSample)
  {
    smart = InputOneSample(t(y[i,]),param,smart,const,initial)
    const[1] = const[1]+1
  }

  smart

}
