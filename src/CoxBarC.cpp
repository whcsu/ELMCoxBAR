#include <vector>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <iostream>
#include <RcppArmadillo.h>
#include "global.h"
#include <RcppNumerical.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Numer;

// Compute the likelihood for the CoxBAR model: 
// y:      vector of responses, where y_i is in [1,...,J_i]
// X:      an R data matrix object, where each element is a matrix of J_i rows and P columns
// params: vector of length P
//
// Returns: a double representing the loglikelihood of the data given the parameters

// [[Rcpp::export]]
double cox_llk_cpp(Rcpp::IntegerVector y,SEXP XData, Rcpp::NumericVector beta,double lambda,Rcpp::NumericVector weight) {

 //convert a dataframe to matrix
	Rcpp:: NumericMatrix X=testDFtoNM(XData);
   //N <- nrow(data) 
	int n = X.nrow();
	arma::mat xbeta=mv_mult(X,beta); 
	arma::colvec expxcol= exp(xbeta.col(0));

    double temp=0; 
	double ll=0; 

  for (int i = 0; i < y.size(); i++) {
	  temp=0;	
    for(int j=i;j<n;j++)
		 temp=temp+expxcol(j);
     ll = ll + y[i]*(xbeta(i,0) - log(temp));   	

  }

  // BAR Regualzied Version
	//elementwise multiplication 
	NumericVector temp1=beta * beta;
	NumericVector temp2=weight * weight+1e-9;
	temp1=temp1/temp2;
	
	//ll-= std::accumulate(temp1.begin(),temp1.end(),0)*0.5*lambda;
	ll-= sum(temp1)*0.5*lambda;
  return ll;
}




// Compute the gradient for the Cox model.  
// Arguments are the same as for llk().  See notes for a derivation.
//
// Returns: an arma::colvec object of length P representing the gradient.
// [[Rcpp::export]]
arma::colvec cox_grad_cpp(Rcpp::IntegerVector y,SEXP XData, Rcpp::NumericVector beta,double lambda,Rcpp::NumericVector weight) {

 //convert a dataframe to matrix
	Rcpp:: NumericMatrix X=testDFtoNM(XData);
   //N <- nrow(data) 
	int n = X.nrow();
	arma::mat xbeta=mv_mult(X,beta); 
	arma::colvec expxcol= exp(xbeta.col(0));

  // Initialize gradient vector
  int P = beta.size(); 
  arma::colvec grad=zeros(P);
  arma::colvec temp1=zeros(P);
  double temp=0;

  for (int i = 0; i < y.size(); i++) {
	  temp=0;
	  temp1=temp1.zeros();   
   for(int j=i;j<n;j++){     

		 temp=temp+expxcol(j);
		 arma::colvec currow=X(j,_);		 
		 temp1=temp1+currow*expxcol(j);
	  }
     arma::colvec currowi=X(i,_);
     grad = grad + y[i]*(currowi - temp1/temp);   	

  }
  /*for( int j=0;j<P;j++)
       Rprintf("%d\n",grad[j]);*/
   arma::colvec temp3 =beta/(weight*weight+1e-9);
  grad = grad-lambda*(temp3); 
  return grad;
}


// Compute the second partial derivatives(Hessian matrix) for the Cox model.  
// Arguments are the same as for llk(). 
// Returns: an arma::mat object of  P*P representing the Hessian matrix.
// [[Rcpp::export]]
arma::mat cox_dgrad_cpp(Rcpp::IntegerVector y,SEXP XData, Rcpp::NumericVector beta,double lambda,Rcpp::NumericVector weight) {

 //convert a dataframe to matrix
	Rcpp:: NumericMatrix X=testDFtoNM(XData);
   //N <- nrow(data) 
	int n = X.nrow();
	int P = beta.size();
	arma::mat xbeta=mv_mult(X,beta); 
	arma::colvec expxcol= exp(xbeta.col(0));

  // Initialize Hessian matrix vector
  arma::mat dgrad(P,P);
  dgrad=dgrad.zeros();
  double temp=0;
  arma::colvec temp1=beta*0;
  arma::mat temp2(P,P);

  for (int i = 0; i<n; i++) {
	  //check user interrupt just in case 
        Rcpp::checkUserInterrupt();
	  temp=0;
	  temp1=temp1.zeros();
	  temp2=temp2.zeros();
	  for(int j=i;j<n;j++){     

		 temp=temp+expxcol(j);
		 arma::colvec currow=X(j,_);
		 
		 temp1=temp1+currow*expxcol(j);
	
		 temp2=temp2+(currow*currow.t())*expxcol(j);
	  }
		dgrad=dgrad+y[i]*(temp2/temp-(temp1*temp1.t())/(temp*temp));
       
   }
 arma::colvec temp3 =lambda/(weight*weight+1e-9);
   dgrad.diag()=dgrad.diag()+temp3;   
return dgrad;

}
