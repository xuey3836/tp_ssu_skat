#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector cut(NumericVector x, NumericVector cutpoint){
  int n = x.size();
  NumericVector g(n);
  for (int i=0;i<n;i++){
    if (x(i)<cutpoint(0)){
      g(i)=0;
    }else if(x(i)<cutpoint(1)){
      g(i)=1;
    }else{
      g(i)=2;
    }
  }
  return g;
}