#include <Rcpp.h>
using namespace Rcpp;

//' Fisher Information for 1-3PL IRT models
//'
//' This function returns fisher information at for a set of items at a particular ability level.
//' @param thetas A vector of values of theta.
//' @param a A vector of a parameters
//' @param b A vector of b parameters
//' @param c A vector of c parameters
//' @export
//' @examples
//' FisherInfo()
// [[Rcpp::export]]





NumericVector FisherInfo(NumericVector thetas, NumericVector a, NumericVector b, NumericVector c) {
  int m = a.size(), N = thetas.size();
  NumericMatrix P(N,m);
  
  for(int j = 0; j < m; j++){
    for(int i = 0; i < N; i++){
      double logit = a[j]*(thetas[i] - b[j]);
      P(i,j) = c[j] + (1-c[j])/(1 + exp(-logit));
    }
  }
  NumericMatrix info(N,m);

  for(int j = 0; j < m; j++){
    for(int i = 0; i < N; i++){
      /*      info(i,j) = ((pow(a[j],2))*(1 - P(i,j))*()); */
      info(i,j) = (pow(a[j],2)*(1 - P(i,j))*pow((P(i,j) - c[j]),2))/(P(i,j)*pow((1-c[j]),2));
    }
  }
  return info;
}

