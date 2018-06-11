



rm(list=ls())
library(SKAT)
data(SKAT.example)
names(SKAT.example)
attach(SKAT.example)
obj<-SKAT_Null_Model(y.b ~ X, out_type="D")
a = SKAT(Z/2, obj,kernel = "linear")


##SSU with covariate vector
n = 100
p=10
k=2
Y = sample(c(0,1),n,prob=c(0.6,0.4),replace = T)
G = matrix(sample(c(0,1,2),size = n*p,replace = T,prob = c(0.1,0.5,0.4)),ncol =p )
library(mvtnorm)
Z = rmvnorm(n=n,mean = c(0,0),sigma = matrix(c(1,0.5,0.5,1),ncol=2))

n<-length(Y)
k<-ncol(G)
k2<-ncol(Z)

#######construction of the score vector and its cov matrix:
if (is.null(Z)){
  ## NO nuisance parameters:
  Xg <- G
  Xbar<-apply(Xg, 2, mean)
  Xgb<-Xg
  for(i in 1:nrow(Xg))
    Xgb[i,]<-Xg[i,]-Xbar
  ##score vector:
  U<-t(Xg) %*% (Y-mean(Y))
  ##cov of the score stats:
  CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
  
} else {
  ## with nuisance parameters:
  tdat1<-data.frame(trait=Y, Z)
  fit1<-glm(trait~.,family="binomial",data=tdat1)
  pis<-fitted.values(fit1)
  Us<-matrix(0, nrow=n, ncol=k)
  for(i in 1:k){
    tdat2<-data.frame(X1=X[,i], Z)
    fit2<-glm(X1~.,data=tdat2)
    X1mus<-fitted.values(fit2)
    Us[, i]<-(Y - pis)*(X[,i] - X1mus)
  }
  U<-apply(Us, 2, sum)
  CovS<-matrix(0, nrow=k, ncol=k)
  for(i in 1:n)
    CovS<-CovS + Us[i,] %*% t(Us[i,])
}


SumSqU<-function(U, CovS){
  if (is.null(dim(CovS))) {# only one-dim:
    Tscore<- sum(U^2 /CovS)
    if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
    pTg1<-as.numeric(1-pchisq(Tscore, 1))
  }
  else {
    #it's possible U=0 and Cov(U)=0:
    if (all(abs(U)<1e-20)) pTg1<-1 else{
      Tg1<- t(U) %*% U
      ##distr of Tg1 is sum of cr Chisq_1:
      cr<-eigen(CovS, only.values=TRUE)$values
      ##approximate the distri by alpha Chisq_d + beta:
      alpha1<-sum(cr*cr*cr)/sum(cr*cr)
      beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
      d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
      alpha1<-as.real(alpha1)
      beta1<-as.real(beta1)
      d1<-as.real(d1)
      pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
    }
  }
  return(pTg1)
}

SumSqUw<-function(U, CovS){
  
  if (is.null(dim(CovS))) {# only one-dim:
    Tg2<- sum(U^2 /CovS)
    if (is.na(Tg2) || is.infinite(Tg2) || is.nan(Tg2)) Tg2<-0
    pTg2<-as.numeric(1-pchisq(Tg2, 1))
  }
  else {
    #it's possible U=0 and Cov(U)=0:
    if (all(abs(U)<1e-20)) pTg2<-1 else{
      diagCovS<-diag(CovS)
      diagCovS<-ifelse(diagCovS>1e-10, diagCovS, 1e-10)
      #Tg2<- t(U) %*%  diag(1/diagCovS) %*% U
      Tg2<- sum(U^2 /diagCovS)
      ##distr of Tg1 is sum of cr Chisq_1:
      cr<-eigen(CovS %*% diag(1/diagCovS), only.values=TRUE)$values
      ##approximate the distri by alpha Chisq_d + beta:
      alpha2<-sum(cr*cr*cr)/sum(cr*cr)
      beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
      d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
      alpha2<-as.real(alpha2)
      beta2<-as.real(beta2)
      d2<-as.real(d2)
      pTg2<-as.numeric(1-pchisq((Tg2-beta2)/alpha2, d2))
    }
  }
  pTg2
}



