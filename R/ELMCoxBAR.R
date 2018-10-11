##'Extreme Learning Machine Cox Model for High Dimensional Survival Analysis
##' @title ELMCoxBAR
##' @param x  The covariates(predictor variables) of training data.
##' @param y  Survival time and censored status of training data. Must be a Surv  \code{survival} object
##' @param Kernel_type Type of kernel matrix. Currently four options avaibable. "RBF_kernel",a RBF kernel;"lin_kernel" , a linear kernel;poly_kernel ,a polynomial kernel;sigmoid_kernel, a sigmoid kernel. Default is "lin_kernel".
##' @param Kernel_para Parameters for different types of kernels. A single value for RBF and linear kernels. A vector for polynomial and sigmoid kernels and progam stops if only a single value is supplied. However, if the vector of values is supplied in the cases of RBF and liner kernels, only the first value will be used. Default is a vector value "c(2,1)".
##' @param penality Currently, penality is defaulted to 0 to train an ELMCoxBAR model. 
##' @param maxiter Maximum values of iterations to update the CoxBAR estimator. Default is 5.
##' @param ... Additional arguments for  glmnet.
##' @return Object of class \code{ELMCoxBAR} with elements
##'   \tabular{ll}{
##'       \code{elmcox}    \tab  A glmnet type model. See \code{glmnet} for details. \cr
##'       \code{trainx} \tab  Training data covariates. \cr
##'          \code{kerneltype} \tab  Type of kernel matrix used in training. kerneltype=1,a RBF kernel;kerneltype=2 , a linear kernel;kerneltype=3 ,a polynomial kernel;kerneltype=4, a sigmoid kernel. \cr
 ##'   \code{Kernel_para} \tab  Parameters used in training. A single value for kerneltype=1 or 2. A vector for kerneltype=3 or 4. \cr
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
##' # Default with lin_kernel 
##' elmsurvmodel=ELMCoxBAR(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]))
##' # with  the RBF kernel
##' elmsurvmodel=ELMCoxBAR(x=trset[,-rii],y=Surv(trset[,rii[1]], 
##' trset[,rii[2]]),Kernel_type="RBF_kernel",Kernel_para=c(2,1))
##' #The predicted linear predictor
##' testprelin=predict(elmsurvmodel,teset[,-c(rii)])
##' @export
ELMCoxBAR <- function(x,y, Kernel_type="lin_kernel", Kernel_para=c(2,1),penality=0, maxiter=5, ...) {

  if (!inherits(y, "Surv"))
    stop("Response must be a 'survival' object - use the 'Surv()' function")

  kplen=length(Kernel_para)

  if(Kernel_type=="RBF_kernel"){
    kerneltype=1
    if(kplen==0||kplen<1){
      stop("Error: Kernel Parameter for RBF_kernel Error!")
    }
  }else if(Kernel_type=="lin_kernel"){
    kerneltype=2
    if(kplen==0||kplen<1){
      stop("Error: Kernel Parameter for lin_kernel Error!")
    }
  }  else if(Kernel_type=="poly_kernel"){
    kerneltype=3
    if(kplen==0||kplen<2){
      stop("Error: Kernel Parameter for poly_kernel Error!")
    }
  }  else if(Kernel_type=="sigmoid_kernel"){
    kerneltype=4
    if(kplen==0||kplen<2){
      stop("Error: Kernel Parameter for sigmoid_kernel Error!")
    }
  }else{
    stop("Error:Unknow kernel types!")
  }

  H = kernmat(x,kerneltype, Kernel_para,NULL)
   
   mylamb=log(nrow(H))
   
   glmcox = glmnet(as.matrix(H), as.matrix(y),family = "cox", lambda=mylamb,alpha = penality)  
   inibeta=as.vector(coef(glmcox))

   fit <- list()
   fit$elmcox=CoxBAR(as.matrix(H), y,weight=inibeta,maxiter)
   fit$trainx = x
   fit$kerneltype = kerneltype
   fit$Kernel_para = Kernel_para

   class(fit)="ELMCoxBAR"
   fit
}





