#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function doubleo an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
CharacterVector defmodel(NumericVector y, NumericVector x){
  double n = y.size();
  double r=0;
  double n1=0;
  double r1=0;
  double s1=0;
  double n2=0;
  double r2=0;
  double s2=0;
  double p1=0;
  double p2=0;
  double delta_p,t2,Z;
  double index;
  for(double i=0;i<n;i++){
    if(y(i)==1){
      r=r+1;
      if(x(i)==1){
        r1=r1+1;
      }else if(x(i)==2){
        r2 =r2+1;
      }
    }else{
      if(x(i)==1){
        s1=s1+1;
      }else if(x(i)==2){
        s2 =s2+1;
      }
    }
  }
  p1 =r1/r;
  p2= r2/r;
  n1 = r1+s1;
  n2= r2+s2;
  delta_p=p2-(p2+p1/2)*(p2+p1/2);
  t2 = (1-n2/n-n1/(2*n))*(n2/n+n1/(2*n));
  Z= sqrt(r)*(delta_p)/t2;
  if(Z> 1.65){
    return "R";
  }else if(Z<-1.65){
    return "D";
  }else{
    return "A";
  }  return Z;
}
  
