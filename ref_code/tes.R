IBS_R = function(X){
  n = nrow(X)
  m= ncol(X)
  K =matrix(1,n,n)
  for (i in 1:n){
    for(j in 1:n){
      temp=0;
      for(l in 1:n){
        diff = 1-abs(X[i,l]-X[j,l])
        temp = temp+ diff
      }
      K[i,j]=K[j,i]= temp/(m)
    }
  }
  return(K)
}

kernel_IBS <-
  function(Z, n, p)
  {
    K = diag(1, n, n)	
    aux = .C("kernel_IBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)))[[4]]
    matrix(aux, nrow=n)
  }
X = matrix(sample(c(0,1,2),9,replace = T),3,3)
X2=X/2
library(Rcpp)
cppFunction("double abss(double x){
            if(x>0){
return x;
            }else{
return -x;
            }
}")

sourceCpp("function/IBS_kernel.cpp")
IBS_kernel(X/2,3,3)
IBS_R(X/2)
kernel_IBS(X,3,3)

library(rbenchmark)
benchmark(c = kernel_IBS(X,n=1000,p=10),
          cpp = IBS_kernel(X=X,1000,10),replications = 1000)[, 1:4]
library(RcppArmadillo)
sourceCpp(file = "mulmat.cpp")
A =matrix(sample(1:100000,10000),1000,10)
B=matrix(sample(1:100000,10000),10,1000)
benchmark(r =t(U) %*%  diag(1/diag(CovS)) %*% U,
          cpp = mul(mul(t(U), diag(1/diag(CovS))),U),replications = 100)[, 1:4]

Rprof()
SNP =gen.MultiSNP(MAF.vec = MAF.vec,sig.mat = sig.mat,n.case=n.case,
                  beta0 = beta0,n.control = n.control,N=N,bet = bet,
                  gen.model=gen.model)
Rprof(NULL)
summaryRprof()
sourceCpp("defmodel.cpp")
def.genmodel(y,x)
defmodel(y,x)
benchmark(r =def.genmodel(y,x),
          cpp = defmodel(y,x),replications = 10000)[, 1:4]
