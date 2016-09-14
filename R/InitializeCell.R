#' InitializeCell
#'
#' initializes parameters in the continuous domain while inputting the first sample
#'
#'
#' @name InitializeCell
#' @param y A one-row matrix, the new sample to be input.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @return The matrix which stores updated parameters over the continuous domain.


InitializeCell <- function (y, param, initial, const)
{
  n = const[2]
  d = const[3]
  K = const[6]

  if (param == TRUE){
    new = matrix(0,nrow=1,ncol=n+2+(1+2*d+2*d^2)*K)
    for (i in 1:K){
      c = initial[1,1+(1+d+d^2)*(i-1)]
      new[1,n+3+(1+2*d+2*d^2)*(i-1)] = c
      mean = initial[1,(2+(1+d+d^2)*(i-1)):(1+d+(1+d+d^2)*(i-1))]
      corr = matrix(initial[1,(2+d+(1+d+d^2)*(i-1)):((1+d+d^2)*i)],nrow=d,ncol=d)
      new[1,(n+4+d+(1+2*d+2*d^2)*(i-1)):(n+3+2*d+(1+2*d+2*d^2)*(i-1))] = mean
      new[1,(n+4+(1+2*d+2*d^2)*(i-1)):(n+3+d+(1+2*d+2*d^2)*(i-1))] = mean*c
      new[1,(n+4+2*d+d^2+(1+2*d+2*d^2)*(i-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(i-1))] = matrix(corr,nrow=1)
      new[1,(n+4+2*d+(1+2*d+2*d^2)*(i-1)):(n+3+2*d+d^2+(1+2*d+2*d^2)*(i-1))] = matrix((corr+t(t(mean))%*%t(mean))*c,nrow=1)
    }
  }
  else{
    new = matrix(0,nrow=1,ncol=n+2+K*d)
    for (i in 1:K){
      mean = initial[1,(1+d*(i-1)):(d*i)]
      new[1,(n+3+d*(i-1)):(n+2+d*i)] = mean
    }
  }

  new[1,1:n] = y[1,1:n]

  new

}
