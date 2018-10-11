# Extreme Learning Machine Cox Model for High Dimensional Survival Analysis
### Introduction
Some interesting recent studies have shown that neural network models are useful alternatives in modelling  survival data when the assumptions of a classical parametric or semiparametric survival model such as the Cox (1972) model are seriously violated. However,  to the best of our knowledge, the plausibility of adapting the emerging extreme learning machine algorithm for single hidden layer feedforward neural networks to survival analysis has not been explored. In this paper, we present a kernel extreme learning machine Cox model regularized by an $L_0$-based broken adaptive ridge (BAR) penalization method. The we demonstrate that the resulting method, referred to as ELMCoxBAR, can outperform some other state-of-art survival prediction methods such as $L_1$- or $L_2$-regularized Cox regression, random survival forest with various splitting rules, and boosted Cox model, in terms of its predictive performance using both simulated and real world datasets. In addition to its good predictive performance,  we illustrate that the proposed method has a key computational advantage over  the above competing methods in terms of computation time efficiency using an  a real-world ultra-high dimensional survival data.

### Installation 
R version >= 3.1 and the latest new Rtools toolchain need to be installed to compile the package. With the "devtools" package, it is easy to install the latest SurvELM R package from Github:
```R
library(devtools)
install_github("whcsu/ELMCoxBAR")
```
#### Sample R code

```R
set.seed(123)
require(ELMCox)
require(survival)
#Lung DATA
data(lung)
lung=na.omit(lung)
lung[,3]=lung[,3]-1
n=dim(lung)[1]
L=sample(1:n,ceiling(n*0.5))
trset<-lung[L,]
teset<-lung[-L,]
rii=c(2,3)
# Default with lin_kernel
elmsurvmodel=ELMCoxBAR(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]))
# with  the RBF kernel
elmsurvmodel=ELMCoxBAR(x=trset[,-rii],y=Surv(trset[,rii[1]],
trset[,rii[2]]),Kernel_type="RBF_kernel",Kernel_para=c(2,1))
#The predicted linear predictor
testprelin=predict(elmsurvmodel,teset[,-c(rii)])
```

### References

Hong Wang and Gang Li. SExtreme Learning Machine Cox Model for High Dimensional Survival Analysis, submitted to Statistics in Medicin
