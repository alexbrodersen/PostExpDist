#include <Rcpp.h>
using namespace Rcpp;

//' Item Characteristic function for 1-3PLM
//'
//' This function returns the probability of getting an item correct under the 1-3PLM.
//' @param thetas A vector of values of theta.
//' @param a A vector of a parameters
//' @param b A vector of b parameters
//' @param c A vector of c parameters
//' @export
//' @examples
//' Probs3plm()
// [[Rcpp::export]]


NumericVector Probs3plm(NumericVector thetas, NumericVector a, NumericVector b, NumericVector c) {
  int m = a.size(), N = thetas.size();
  NumericMatrix P(N,m);
  
  for(int j = 0; j < m; j++){
    for(int i = 0; i < N; i++){
      double logit = a[j]*(thetas[i] - b[j]);
      P(i,j) = c[j] + (1-c[j])/(1 + exp(-logit));
    }
  }
  return P;  
}

