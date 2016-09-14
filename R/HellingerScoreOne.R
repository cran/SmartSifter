#' HellingerScoreOne
#'
#' calculates the Hellinger score after inputting one sample y
#'
#'
#' @name HellingerScoreOne
#' @param y A one-row matrix, the new sample to be input.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @return The Hellinger score after inputting one sample y.



HellingerScoreOne <- function (y, param = TRUE, smart, const, initial)
{
  r = const[5]
  K = const[6]
  n = const[2]
  d = const[3]

  newsmart = InputSample (y, param, smart, const, initial)
  #you get the updated matrix of parameters

  oldsmart = smart

  if (dim(newsmart)[1] > dim(smart)[1]){
    new = InitializeCell(y, param, initial, const)
    oldsmart = rbind(smart,new)
  }

  numOfCell = dim(newsmart)[1]

  term1 = 0
  term2 = 0

  #in parametric version, correlation matrix is stored in column-major
  if (param == TRUE){
    for (i in 1:numOfCell){
      tmp = sqrt(oldsmart[i,n+2]*newsmart[i,n+2])
      term1 = term1 + tmp
      dh = 0
      for (j in 1:K){
        corrold = matrix(oldsmart[i,(n+4+2*d+d^2+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(j-1))], nrow=d,ncol=d)
        corrnew = matrix(newsmart[i,(n+4+2*d+d^2+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(j-1))], nrow=d,ncol=d)
        meanold = matrix(oldsmart[i,(n+4+d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+(1+2*d+2*d^2)*(j-1))],ncol=1)
        meannew = matrix(newsmart[i,(n+4+d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+(1+2*d+2*d^2)*(j-1))],ncol=1)
        I = diag(1,nrow=d,ncol=d)
        invcorrold = solve(corrold,I)
        invcorrnew = solve(corrnew,I)
        tmp1 = 2*(det((invcorrold+invcorrnew)/2))^(-1/2)
        tmp1 = tmp1/((det(corrold))^(1/4)*(det(corrnew))^(1/4))
        tmp2 = 1/2*t(invcorrold%*%meanold+invcorrnew%*%meannew)%*%solve(invcorrold+invcorrnew,I)%*%(invcorrold%*%meanold+invcorrnew%*%meannew)
        tmp2 = exp(tmp2)
        tmp3 = 1/2*(t(meanold)%*%invcorrold%*%meanold+t(meannew)%*%invcorrnew%*%meannew)
        tmp3 = exp(tmp3)
        appdh = 2-tmp1*tmp2*tmp3
        cold = oldsmart[i,n+3+(1+2*d+2*d^2)*(j-1)]
        cnew = newsmart[i,n+3+(1+2*d+2*d^2)*(j-1)]
        dh = dh + (sqrt(cold)-sqrt(cnew))^2
        dh = dh + (cold+cnew)/2*appdh
      }
      term2 = term2 + tmp*dh
    }
    score = 1/r^2*(2-2*term1+term2)
  }

  else{
    sigma2 = const[7]
    for (i in 1:numOfCell){
      tmp = sqrt(oldsmart[i,n+2]*newsmart[i,n+2])
      term1 = term1 + tmp
      dh = 0
      for (j in 1:K){
        corrold = diag(sigma2,nrow=d,ncol=d)
        corrnew = diag(sigma2,nrow=d,ncol=2)
        meanold = matrix(oldsmart[i,(n+3+(j-1)*d):(n+2+j*d)],ncol=1)
        meannew = matrix(newsmart[i,(n+3+(j-1)*d):(n+2+j*d)],ncol=1)
        I = diag(1,nrow=d,ncol=d)
        invcorrold = solve(corrold,I)
        invcorrnew = solve(corrnew,I)
        tmp1 = 2*(det((invcorrold+invcorrnew)/2))^(-1/2)
        tmp1 = tmp1/((det(corrold))^(1/4)*(det(corrnew))^(1/4))
        tmp2 = 1/2*t(invcorrold%*%meanold+invcorrnew%*%meannew)%*%solve(invcorrold+invcorrnew,I)%*%(invcorrold%*%meanold+invcorrnew%*%meannew)
        tmp2 = exp(tmp2)
        tmp3 = 1/2*(t(meanold)%*%invcorrold%*%meanold+t(meannew)%*%invcorrnew%*%meannew)
        tmp3 = exp(tmp3)
        appdh = 2-tmp1*tmp2*tmp3
        cold = 1/K
        cnew = 1/K
        dh = dh + (sqrt(cold)-sqrt(cnew))^2
        dh = dh + (cold+cnew)/2*appdh
      }
      term2 = term2 + tmp*dh
    }
    score = 1/r^2*(2-2*term1+term2)
  }

  score

}
