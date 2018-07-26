#include <Rcpp.h>
#include <cmath>        // std::abs
#include <math.h>       /* isinf, sqrt */
using namespace Rcpp;


//' Maximum Likelihood Scoring for Item Response Theory
//'
//' This function estimates IRT ability parameters using Fisher Scoring and then interval bisection if that fails.
//' @param irvs A matrix of responses with persons on rows and items on columns. No missing data allowed.
//' @param initThetas A vector of initial values of theta.
//' @param a A vector of a parameters
//' @param b A vector of b parameters
//' @param c A vector of c parameters
//' @param maxit The maximum number of iterations for Fisher Scoring.
//' @param crit A value such that the difference between successive thetea iterations indicates convergence.
//' @export
//' @examples
//' MLscoring()
// [[Rcpp::export]]



List MLscoring(IntegerMatrix irvs, NumericVector initThetas, NumericVector a, NumericVector b, NumericVector c,int maxit = 25, double crit = 1e-4) {
  
  int m = a.size(), N = initThetas.size();  
  NumericVector P(m);
  NumericVector thetas(N);
  for(int i = 0; i < N; i++){
    thetas[i] = NA_REAL;
  }

  NumericVector se(N);
  
  
  for(int i = 0; i < N; i++){
    int total = 0;

    for (int j = 0; j < m; j++) {
      total += irvs(i, j);
    }
    
    if(total == 0){
      continue;
    } else if(total == m){
      continue;
    }

    double thetaTemp = initThetas[i];
    double deltaN = 0;
    double deltaD = 0;
    for(int r = 0; r < maxit; r++){
      deltaN = 0;
      deltaD = 0;
    double delta;      
      for(int j = 0; j < m; j++){
	  double logit = a[j]*(thetaTemp - b[j]);
	  P[j] = c[j] + (1-c[j])/(1 + exp(-logit));
      }
      
      for(int j = 0; j < m; j++){
	deltaN += (a[j]*(irvs(i,j) - P[j])*(P[j] - c[j]))/(P[j]*(1-c[j]));
	deltaD += ((a[j]*a[j]*(1-P[j]))*pow((P[j] - c[j]),2))/(P[j]*pow((1-c[j]),2));
	  }

      delta = deltaN/deltaD;
      if(std::abs(delta) > 2){
	break;
      } else if(std::isinf(delta) || std::isnan(delta)){
	break;
      }

      thetaTemp = thetaTemp + delta ;
      if(std::abs(delta) < crit){
	thetas[i] = thetaTemp;
	break;
      }
    }
    
    if(NumericVector::is_na(thetas[i])){
      double L = -5;
      double R = 5;
      double M = (R + L)/2; 
      int iters = ceil(log2((R - L)/crit) - 1);      
      NumericVector LMR =  NumericVector::create(L, M, R);      
      int K = LMR.size();
      NumericMatrix tP(K,m);
 
      for(int r = 0; r < iters; r ++){
	M = (R + L)/2;
	LMR =  NumericVector::create(L, M, R);      
	for(int j = 0; j < m; j++){
	  for(int k = 0; k < K; k++){
	    double logit = a[j]*(LMR[k] - b[j]);
	    tP(k,j) = c[j] + (1-c[j])/(1 + exp(-logit));
	  }
	}
	NumericVector dlnL = NumericVector::create(0,0,0);
	for(int k = 0; k < K; k++){
	  for(int j = 0; j < m; j++){
	    dlnL[k] += a[j]*(irvs(i,j) - tP(k,j))*(tP(k,j) - c[j])/(tP(k,j)*(1-c[j]));
	  }
	}
	  if(dlnL[0]*dlnL[2] > 0){
	    break;
	  }
	  
	  IntegerVector s(K);
	  for(int k = 0; k < K; k ++){
	    if(dlnL[k] > 0){
	      s[k] = 1;
	    } else if(dlnL[k] < 0){
	      s[k] = 0;
	    }
	  }	  
	  if(s[0] + s[1] == 1){
	    R = M;
	  } else if(s[0] + s[1] != 1){
	    L = M;
	  }
	  if(r == (iters - 1)){
	    thetas[i] = (L + R)/2;
	  }
	  
      }
    }


    NumericVector sP(m);
    for(int j = 0; j < m; j++){
      double logit = a[j]*(thetas[i] - b[j]);
      sP(j) = c[j] + (1-c[j])/(1 + exp(-logit));
    }

    double info = 0;
    
    for(int j = 0; j < m; j++){
      info += (pow(a[j],2)*(1 - sP(j))*pow((sP(j) - c[j]),2))/(sP(j)*pow((1-c[j]),2));
    }
    se(i) = sqrt(1/info);
  }
  
  return List::create(Rcpp::Named("thetas") = thetas,Rcpp::Named("se") = se);
}


