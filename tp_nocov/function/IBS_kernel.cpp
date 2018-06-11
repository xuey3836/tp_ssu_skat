#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericMatrix IBS_kernel(NumericMatrix X, int n, int m) {
  double diff;
  NumericMatrix K(n,n);
  K.fill_diag(1);
  double temp;
  for (int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      temp=0;
      for(int l=0; l<m; l++){
        diff = 1-std::abs(X(i,l)-X(j,l));
        temp = temp+ diff;
      }
      K(i,j)=K(j,i)= temp/m;
    }
  }
  return K;
}