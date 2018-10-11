##'Cox's Regression with broken adaptive ridge (CoxBAR)
##' @title CoxBAR
##' @param trainx  The covariates(predictor variables) of training data.
##' @param y  Survival time and censored status of training data. Must be a Surv  \code{survival} object
##' @param weight In ELMCoxBAR, we set this to a random Cox-Lasso estimate. 
##' @param maxiter Maximum values of iterations to update the CoxBAR estimator. Default is 5.
##' @param standardize Logical flag for trainx variable standardization, prior to fitting the model sequence. Default is standardize=TRUE
##' @return Object of class \code{CoxBAR} with elements
##'   \tabular{ll}{
##'       \code{meanx} \tab  Mean values of original trainx if standardization is TRUE. \cr
##'       \code{sdx} \tab  Standard deviation values of original trainx if standardization is TRUE. \cr
##'       \code{standardize} \tab  The standardization status. \cr
##'       \code{beta}    \tab   The point  estimates of \eqn{\beta}. \cr
##'       \code{logLik} \tab  Log Likelihood. \cr
##'        
##'   }
##' @author Hong Wang
##' @examples
##' set.seed(123)
##' require(ELMCox)
##' require(survival)
##' #Lung DATA
##' data(lung)
##' lung=na.omit(lung)
##' lung[,3]=lung[,3]-1
##' n=dim(lung)[1]
##' L=sample(1:n,ceiling(n*0.5))
##' trset<-lung[L,]
##' teset<-lung[-L,]
##' rii=c(2,3)
##' # A randon weight for illustration purpose.
##' p=dim(lung)[2]-2
##' myweight=rep(0.5,p)
##' coxbarmodel=CoxBAR(trainx=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]),weight=myweight)
##' @export
CoxBAR <- function(trainx,y,weight=NULL,maxiter=5,standardize=TRUE) {
  if (!inherits(y, "Surv"))
    stop("Response must be a 'survival' object - use the 'Surv()' function")
  object <- list()
#message(maxiter)
if (standardize) {
    sdx <- apply(trainx,2,sd)
    sdx <- ifelse(sdx == 0,1,sdx)
    meanx <- apply(trainx,2,mean)
    trainx <- scale(trainx,center=meanx,scale=sdx)
    object$meanx <- meanx
    object$sdx <- sdx
    object$standardize <- TRUE
  } else {
    object$standardize <- FALSE
  }
  survtime = y[, 1]
  status <- y[, 2]
  ot <- order(survtime)
  status <- status[ot]
  survtime <- survtime[ot]
  trainx<-trainx[ot,]
  p=ncol(trainx)
  if(is.null(weight)){
     weight=rep(0,p)
  }
  #nlm(f=CoxLLKC,g=CoxGrad,p=beta,trainx,status,mylambda)

  #weight using ridge regression
  #initial values for coefficient
  beta=weight
  lambda=log(length(status))
  iter=0
  while (iter<maxiter){
   #betapre=coef(optimx(beta,fn=CoxLLKC,gr=CoxGradC,hess=CoxDGradC,trainx=trainx,status=status,mylambda=lambda,weight=weight,method="Nelder-Mead"))
 #  betapre=coef(optimx(beta,fn=CoxLLKC,gr=CoxGradC,trainx=trainx,status=status,mylambda=lambda,weight=weight,method="Nelder-Mead"))
   betares=optim(beta,fn=CoxLLKC,trainx=trainx,status=status,mylambda=lambda,weight=weight)
   betapre=betares$par
   #betapre=optim(beta,fn=CoxLLKC,trainx=trainx,status=status,mylambda=lambda,weight=weight)$par
  dx = max(abs(betapre-beta))
  #change from 1e-5 to 1e-9
  if(dx <= 1.0e-5) {iter = maxiter}
   else  {iter=iter+1}
  weight=beta
  beta=betapre
  #message(paste("iter=,beta1=",iter,beta[1]))
  #betapre=optim(beta,fn=CoxLLKC,gr=CoxGrad,trainx=trainx,status=status,mylambda=lambda,method="BFGS")

 }
 object$beta=beta
 object$logLik=betares$value
 object
}


CoxLLK <- function(pv,trainx,status,mylambda) {
  a <- trainx %*% as.matrix(pv)
  b <- log(rev(cumsum(rev(exp(a)))))
  -sum((a - b)[status==1])
}
##' @export
CoxLLKC <- function(pv,trainx,status,mylambda,weight) {
  -cox_llk_cpp(status,trainx,pv,mylambda,weight)
}


##' @export
CoxGradC <- function(pv,trainx,status,mylambda,weight) {
  -cox_grad_cpp(status,trainx,pv,mylambda,weight)
}

##' @export
CoxDGradC <- function(pv,trainx,status,mylambda,weight) {
  cox_dgrad_cpp(status,trainx,pv,mylambda,weight)
}
