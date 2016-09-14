#' LogLossOne
#'
#' calculates the logarithmic loss after inputting one sample
#'
#'
#' @name LogLossOne
#' @param y A one-row matrix, the new sample to be input.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @return The logarithmic loss after inputting one sample y.
#' @importFrom mvtnorm pmvnorm




LogLossOne <- function (y, param, smart, const, initial)
{

  n = const[2]
  d = const[3]
  K = const[6]

  numOfCell = dim(smart)[1]

  index = WhichCell (y, param, smart, const)

  x = y[1,(n+1):(n+d)]

  if (index > numOfCell) { result = Inf }
  else {
    term1 = smart[index,n+2]
    term2 = 0
    if (param == TRUE){
      for (i in 1:K){
        tmp = pmvnorm(as.vector(t(x)), mean=as.vector(matrix(smart[index,(n+4+d+(1+2*d+2*d^2)*(i-1)):(n+3+2*d+(1+2*d+2*d^2)*(i-1))],ncol=1)), sigma=matrix(smart[index,(n+4+2*d+d^2+(1+2*d+2*d^2)*(i-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(i-1))],nrow=d,ncol=d))
        tmp = tmp[1]
        tmp = tmp * smart[index,n+3+(1+2*d+2*d^2)*(i-1)]
        term2 = term2 + tmp
      }
    }
    else {
      sigma2 = const[7]
      for (i in 1:K){
        tmp = pmvnorm(as.vector(t(x)), mean=as.vector(matrix(smart[index,(n+3+d*(i-1)):(n+2+d*i)],ncol=1)), sigma=diag(sigma2,nrow=d,ncol=d))
        tmp = tmp[1]
        tmp = tmp/K
        term2 = term2 + tmp
      }
    }
    result = -log(term1*term2)
  }

  result

}
