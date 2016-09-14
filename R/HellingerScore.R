#' HellingerScore
#'
#' calculates the Hellinger score after inputting sample y (can be more than one)
#'
#' @name HellingerScore
#' @param y A matrix, the new sample to be input.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @return The Hellinger score after inputting sample y.





HellingerScore <- function (y, param = TRUE, smart, const, initial)
{
  NumOfSample = dim(y)[1]

  score = matrix(0,nrow=1,ncol=NumOfSample)

  for (i in 1:NumOfSample){
    score[1,i] = HellingerScoreOne (t(y[i,]),param,smart,const,initial)
    smart = InputSample (t(y[i,]),param,smart,const,initial)
    const[1] = const[1]+1
  }

  score

}

