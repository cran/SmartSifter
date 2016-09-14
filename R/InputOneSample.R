#' InputOneSample
#'
#' updates parameters after inputting one sample
#'
#'
#' @name InputOneSample
#' @param y A one-row matrix, the new sample to be input.
#' @param param A logical scalar. If TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @return The matrix which stores updated parameters over the continuous domain.
#' @importFrom mvtnorm pmvnorm
#' @importFrom rootSolve multiroot


InputOneSample <- function (y, param = TRUE, smart, const, initial)
{
  N = const[1] + 1
  n = const[2]
  d = const[3]
  rh = const[4]
  r = const[5]
  K = const[6]
  sigma2 = const[7] #for parametric version
  alpha = const[7] #for non-parametric version

  numOfCell = dim(smart)[1]

  index = WhichCell(y, param, smart, const) #This function returns the row number to which the new sample belongs

  #if the cell does not exist, create a new one
  if (index == numOfCell + 1){
    new = InitializeCell (y,param,initial,const)
    smart = rbind(smart,new)
  }

  numOfCell = dim(smart)[1]

  #update the discrete-property-related values
  for (i in 1:numOfCell){
    if (i == index){
      smart[i,n+1] = (1 - rh) * smart[i,n+1] + 1
      #the (n+1)th column stores T
      smart[i,n+2] = smart[i,n+1]*rh/(1-(1-rh)^N)
      #the (n+2)th column stores p
    }
    else {
      smart[i,n+1] = (1 - rh) * smart[i,n+1]
      smart[i,n+2] = smart[i,n+1]*rh/(1-(1-rh)^N)
    }
  }

  #update the continuous-property-related values
  x = as.matrix(y[,(n+1):(n+d)])

  #if is parametric version, update parameters of continuous domain in smart[index,]
  #in parametric version, the correlation matrix is stored in column-major
  if (param == TRUE){
    sum = 0
    for (j in 1:K) #for every gaussian function
    {
      cov = matrix(smart[index,(n+4+2*d+d^2+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(j-1))],nrow=d,ncol=d)
      mean = smart[index,(n+4+d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+(1+2*d+2*d^2)*(j-1))]
      numerator = pmvnorm(as.vector(t(x)), mean=mean, sigma=cov)
      numerator = numerator[1]
      numerator = numerator * smart[index,n+3+(1+2*d+2*d^2)*(j-1)]
      sum = sum + numerator
    }
    for (j in 1:K) #for every gaussian function, update parameters
    {
      numerator = pmvnorm(as.vector(t(x)), mean=as.vector(matrix(smart[index,(n+4+d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+(1+2*d+2*d^2)*(j-1))],ncol=1)), sigma=matrix(smart[index,(n+4+2*d+d^2+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(j-1))],nrow=d,ncol=d))
      numerator = numerator[1]
      numerator = numerator * smart[index,n+3+(1+2*d+2*d^2)*(j-1)]
      tmp = (1-alpha*r)*numerator/sum + alpha*r/K
      smart[index,n+3+(1+2*d+2*d^2)*(j-1)] = (1-r)*smart[index,n+3+(1+2*d+2*d^2)*(j-1)] + r*tmp   #update pi
      smart[index,(n+4+(1+2*d+2*d^2)*(j-1)):(n+3+d+(1+2*d+2*d^2)*(j-1))] = (1-r)*smart[index,(n+4+(1+2*d+2*d^2)*(j-1)):(n+3+d+(1+2*d+2*d^2)*(j-1))]+r*tmp*x   #update miu__
      corrba = matrix(smart[index,(n+4+2*d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+d^2+(1+2*d+2*d^2)*(j-1))],nrow=d,ncol=d)
      corrba = (1-r)*corrba + r*tmp*x%*%t(x)
      smart[index,(n+4+2*d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+d^2+(1+2*d+2*d^2)*(j-1))] = matrix(corrba, nrow=1)    #update sigma__
      smart[index,(n+4+d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+(1+2*d+2*d^2)*(j-1))] = smart[index,(n+4+(1+2*d+2*d^2)*(j-1)):(n+3+d+(1+2*d+2*d^2)*(j-1))]/smart[index,n+3+(1+2*d+2*d^2)*(j-1)]   #update miu
      corr = smart[index,(n+4+2*d+d^2+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(j-1))]
      mytmp = smart[index,n+3+(1+2*d+2*d^2)*(j-1)]
      mytmp2 = as.matrix(smart[index,(n+4+d+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+(1+2*d+2*d^2)*(j-1))])
      corr = corrba/mytmp - mytmp2%*%t(mytmp2)
      smart[index,(n+4+2*d+d^2+(1+2*d+2*d^2)*(j-1)):(n+3+2*d+2*d^2+(1+2*d+2*d^2)*(j-1))] = matrix(corr,nrow=1)   #update sigma
    }
  }

  #if is non-parametric version, update parameters of continuous domain in smart[index,]
  else {
    f = matrix(0,nrow=K,ncol=d)
    C = matrix(0,nrow=K,ncol=d)
    B = matrix(0,nrow=K,ncol=d)
    model = function (myq){
      for (k in 1:K){
        for (l in 1:d){
          f[k,l] = 0
          for (j in 1:K){
            for (m in 1:d){
              C[j,m] = delta(m,l)-(smart[index,n+2+(k-1)*d+l]-smart[index,n+2+(j-1)*d+l])*(smart[index,n+2+(k-1)*d+m]-smart[index,n+2+(j-1)*d+m])/(2*sigma2)
              tmp = norm(smart[index,(n+3+(k-1)*d):(n+2+k*d)]-smart[index,(n+3+(j-1)*d):(n+2+j*d)],'2')
              C[j,m] = C[j,m]*exp(-tmp*tmp/(4*sigma2))
              f[k,l] = f[k,l] + C[j,m]*myq[(j-1)*d+m]
            }
          }
          tmp = 0
          for (i in 1:K){
            tmp1 = norm(smart[index,(n+3+(i-1)*d):(n+2+i*d)]-smart[index,(n+3+(k-1)*d):(n+2+k*d)],'2')
            tmp = tmp + (smart[index,n+2+(i-1)*d+l]-smart[index,n+2+(k-1)*d+l])*exp(-tmp1*tmp1/(4*sigma2))
          }
          tmp1 = norm(x-smart[index,(n+3+(k-1)*d):(n+2+k*d)],'2')
          B[k,l] = r * ( K*(x[l,1]-smart[index,n+2+(k-1)*d+l])*exp(-tmp1*tmp1/(4*sigma2)) - tmp )
          f[k,l] = f[k,l]-B[k,l]
        }
      }
      result = matrix(0,nrow=1,ncol=K*d)
      for (k in 1:K){
        for (l in 1:d){
          result[d*(k-1)+l] = f[k,l]
        }
      }
      result
    }
    solution = multiroot(f=model,start=matrix(0,nrow=1,ncol=K*d))
    smart[index,(n+3):(n+2+K*d)] = smart[index,(n+3):(n+2+K*d)] + solution$root
    #finish updating q
  }

  smart

}
