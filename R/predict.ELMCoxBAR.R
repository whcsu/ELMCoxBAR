##' Predicting  with an ELMCoxBAR model
##' @title Predict.ELMCoxBAR
##' @param object  An object that inherits from class ELMBJEN.
##' @param testx  A data frame in which to look for variables with which to predict. 
##' @return produces a vector of predictions or a matrix of predictions
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
##' # Default with lin_kernel 
##' elmsurvmodel=ELMCoxBAR(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]))
##' # with  the RBF kernel
##' elmsurvmodel=ELMCoxBAR(x=trset[,-rii],y=Surv(trset[,rii[1]], 
##' trset[,rii[2]]),Kernel_type="RBF_kernel",Kernel_para=c(2,1))
##' #The predicted linear predictor
##' testprelin=predict(elmsurvmodel,teset[,-c(rii)])
##' @export
predict.ELMCoxBAR <- function(object, testx) {
  Kernel_para = object$Kernel_para
  kerneltype = object$kerneltype
  trainx = object$trainx
  beta=object$elmcox$beta
 # message (paste("My var",object$Kernel_para))
  H = kernmat(trainx,kerneltype, Kernel_para,testx)

  if (object$elmcox$standardize) {
    H=scale(H, center=object$elmcox$meanx,scale=object$elmcox$sdx)
    }
  elmcoxpre = H %*%beta

  return(elmcoxpre)
}





