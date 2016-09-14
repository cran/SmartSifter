#' Test
#'
#' input new samples and compute their related scores (to detect possible outliers)
#'
#' @name Test
#' @param y A matrix, the training set.
#' @param param A logical scalar, if TRUE, the model is in parametric version, otherwise, a non-parametric one.
#' @param smart A matrix, stores all the parameters over the continuous domain.
#' @param const A numeric vector, specifies the value of all global variables, if param = T, then const = c(N,n,d,rh,r,K,alpha); if param=FALSE, then const = c(N,n,d,rh,r,K,sigma_sqare).
#' @param initial A numeric vector, specifies the initial value of parameters over the continuous domain, if param = T, initial = c(pi_1,mean_1,cov_1, ..., pi_K, mean_K,cov_K), if param = F, initial = c(q1,q2, ..., qK).
#' @return List of updated parameters.
#' @examples
#' ## The parametric version
#' initial=matrix(c(0.5,0,0,1,0,0,1,0.5,1,1,1,0,0,1),nrow=1)
#' const = c(0,1,2,0.1,0.1,2,2)
#' param=TRUE
#' y=matrix(c(1,3,1,0,1,1),nrow=2)
#' smart = Train(y,param,const,initial)$smart
#' const[1] = Train(y,param,const,initial)$N
#' y=matrix(c(2,1,0),nrow=1)
#' smart = Test(y,param,smart,const,initial)$smart
#' HellingerScore = Test(y,param,smart,const,initial)$HellingerScore
#' LogLoss = Test(y,param,smart,const,initial)$LogLoss
#' const[1] = Test(y,param,smart,const,initial)$N

#' ##The nonparametric version
#' param=FALSE
#' const = c(0,1,2,0.1,0.1,2,1)
#' initial = matrix(c(0,0,1,1),nrow=1)
#' y=matrix(c(1,3,1,0,1,1),nrow=2)
#' smart = Train(y,param,const,initial)$smart
#' const[1] = Train(y,param,const,initial)$N
#' y=matrix(c(2,1,0),nrow=1)
#' smart = Test(y,param,smart,const,initial)$smart
#' HellingerScore = Test(y,param,smart,const,initial)$HellingerScore
#' LogLoss = Test(y,param,smart,const,initial)$LogLoss
#' const[1] = Test(y,param,smart,const,initial)$N
#' @export



Test <- function (y, param=TRUE, smart, const, initial)
{
  smart = InputSample (y, param, smart, const, initial)
  const = UpdateConst (y, const)

  N = const[1]
  n = const[2]
  d = const[3]
  rh = const[4]
  r = const[5]
  K = const[6]

  if (param == TRUE){
    alpha = const[7]
    score = HellingerScore(y,param,smart,const,initial)
    loss = LogLoss(y,param,smart,const,initial)
    out = list (N=N,n=n,d=d,rh=rh,r=r,K=K,alpha=alpha,smart=smart,HellingerScore=score,LogLoss=loss)
  }
  else{
    sigma2 = const[7]
    score = HellingerScore(y,param,smart,const,initial)
    loss = LogLoss(y,param,smart,const,initial)
    out = list (N=N,n=n,d=d,rh=rh,r=r,K=K,sigma2=sigma2,smart=smart,HellingerScore=score,LogLoss=loss)
  }

  out

}

