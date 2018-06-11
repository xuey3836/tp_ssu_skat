
##compute p-value of tpSSU,tpSKAT, SSU and SKAT
#for example
y = sample(c(0,1),size = 1000,prob = c(0.4,0.6),replace = T)
G = matrix(sample(c(0,1,2),size=10000,prob=c(0.1,0.5,0.4),replace = T),nrow = 1000,ncol=10)
Z = as.matrix(rnorm(1000))
## have covariate
system.time(pvalue_tp(y=y,G=G,Z=Z,B=1000))##i5-6200 CPU @2.30GHz, about 40s
## no convariate
pvalue_tp(y=y,G=G,B=1000)##i5-6200 CPU @2.30GHz, about 17s

library(Rcpp)
library(MASS)
##compute the pvalue of tpssu tpskat, ssu,skat
pvalue_tp<-function(y,G,Z=NULL,B=200){
  #two-SSU
  p = ptwostep(y=y,G=G,Z=Z,B=B)
  # p = pboot(y=y,G=G,Z=Z,B=B)
  
  # SSU and SKAT
  p["SSU"] = SSUw(y=y, G = G,Z=Z)
  p["SKAT"]= SKAT.m(y=y,G=G,Z=Z,perm = 0)$asym.pval
  p
}

#SSU and SKAT statistics
stat_ssu_skat<-function(y,G, Z=NULL){
  n<-length(y)
  m<-ncol(G)
  k<-ncol(Z)
  #######construction of the score vector and its cov matrix:
  if (is.null(Z)){
    ## NO nuisance parameters:
    Gc<-apply(G, 2, function(x)(x-mean(x)) )
    y.new = y - mean(y)
    ##score vector:
    U<-t(G) %*% y.new
    ##cov of the score stats:
    CovS<- mean(y)*(1-mean(y))*(t(Gc) %*% Gc)
  }else {
    tdat1<-data.frame(trait=y, Z)
    fit1<-glm(trait~.,family="binomial",data=tdat1)
    pis<-fitted.values(fit1)
    U <- t(G) %*% (y-pis)
    X = cbind(G,Z)
    Xc <- apply(X, 2, function(x)(x-mean(x)) )
    ##cov of the score stats:
    ######Cov of the WHOLE score stat vector:
    CovS0<-t(Xc) %*% diag(pis*(1-pis)) %*% Xc
    ##cov of the score stats:
    CovS<- ginv(ginv(CovS0)[1:m, 1:m])
    y.new = y - pis
  }
  #weighted one:
  Tg2<- t(U) %*%  diag(1/diag(CovS)) %*% U
  ##skat.stat
  # K1= kernel_IBS(2*X, n, m)
  if (max(G)==2){
    K= IBS_kernel(G/2,n = n,m = m)
  }else{
    K= IBS_kernel(G,n = n,m = m)
  }
  Q = t(y.new) %*% K %*% y.new/2
  skat.stat = as.numeric(Q)
  return(list(SSU= as.numeric(Tg2),SKAT =as.numeric(skat.stat)))
}


TwoStep<-function(y,G,Z=NULL){
  STAT=NULL
  m=ncol(G)
  ##select gentic model
  gen.model=apply(G,2,function(x){
    # def.genmodel(y=y,x=x)
    defmodel(y=y,x=x)
  })
  GX=sapply(1:m,function(j){if(gen.model[j]=="R"){
    as.numeric(G[,j]>1)
  }else if(gen.model[j]=="D"){
    as.numeric(G[,j]>0)
  }else{G[,j]/2}})
  
  ## compute statistics
  stat=stat_ssu_skat(y=y,G=GX,Z=Z)
  STAT["tpSSU"] =stat$SSU
  STAT["tpSKAT"] = stat$SKAT
  return(STAT)
}


# B=1000,permutation
ptwostep<-function(y,G,Z=NULL,B){
  G = as.matrix(G)
  Stat=TwoStep(y = y,G=G,Z=Z)
  nc=length(Stat)
  T2=matrix(0,ncol=nc,nrow = B)
  for(i in 1:B){
    if(i %% 100 ==0){
      cat("the premutation step is",i,"\n")
    }
    perm.sample = sample(1:length(y))
    y.perm = y[perm.sample]
    Z.perm = Z[perm.sample,]
    T2[i,]=TwoStep(y = y.perm,G=G,Z=Z.perm)
  }
  p = sapply(1:nc,function(j){mean(T2[,j]>Stat[j])})
  names(p)<-names(Stat)
  return(p)
}

#HWET to select the genetic model
cppFunction('CharacterVector defmodel(NumericVector y, NumericVector x){
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
            ')

##IBS kernel
cppFunction('NumericMatrix IBS_kernel(NumericMatrix X, int n, int m) {
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
            }')


##
SSUw<-function(y, G,Z=NULL){
  n<-length(y)
  m<-ncol(G)
  k<-ncol(Z)
  #######construction of the score vector and its cov matrix:
  if (is.null(Z)){
    ## NO nuisance parameters:
    Gc<-apply(G, 2, function(x)(x-mean(x)) )
    y.new = y - mean(y)
    ##score vector:
    U<-t(G) %*% y.new
    CovS<- mean(y)*(1-mean(y))*(t(Gc) %*% Gc)
  }else {
    tdat1<-data.frame(trait=y, Z)
    fit1<-glm(trait~.,family="binomial",data=tdat1)
    pis<-fitted.values(fit1)
    U <- t(G) %*% (y-pis)
    X = cbind(G,Z)
    Xc <- apply(X, 2, function(x)(x-mean(x)) )
    ##cov of the score stats:
    ######Cov of the WHOLE score stat vector:
    CovS0<-t(Xc) %*% diag(pis*(1-pis)) %*% Xc
    ##cov of the score stats:
    CovS<- ginv(ginv(CovS0)[1:m, 1:m])
  }
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
      alpha2<-as.numeric(alpha2)
      beta2<-as.numeric(beta2)
      d2<-as.numeric(d2)
      pTg2<-as.numeric(1-pchisq((Tg2-beta2)/alpha2, d2))
    }
  }
  pTg2
}

##SKAT
SKAT.m = function (y, G, Z=NULL,perm =0) 
{
  G = as.matrix(G)
  n = nrow(G)
  m = ncol(G)
  if (max(G)==2){
    K= IBS_kernel(G/2,n = n,m = m)
  }else{
    K= IBS_kernel(G,n = n,m = m)
  }
  if (is.null(Z)){
    mu = rep(mean(y),n)
  }else{
    tdat1<-data.frame(trait=y, Z)
    fit1<-glm(trait~.,family="binomial",data=tdat1)
    mu <-fitted.values(fit1)
  }
  y.new = y - mu
  Q = t(y.new) %*% K %*% y.new/2
  V = diag(mu * (1 - mu))
  X = cbind(rep(1,n),Z)
  P = V - V %*% X %*% ginv(t(X) %*% V %*% X) %*% 
    t(X) %*% V
  PK = P %*% K
  muQ = sum(diag(PK))/2
  itau = sum(PK * t(PK))/2
  itausig = sum(PK * t(P))/2
  isigma = sum(P^2)/2
  iest = itau - ((itausig^2)/isigma)
  skale = iest/(2 * muQ)
  degfr = 2 * (muQ^2)/iest
  skat.stat = as.numeric(Q)
  asym.pval = 1 - pchisq(skat.stat/skale, df = degfr)
  perm.pval = NA
  if (perm > 0) {
    G.perm = rep(0, perm)
    for (i in 1:perm) {
      perm.sample = sample(1:length(y))
      y.perm = y[perm.sample]
      tdat1<-data.frame(trait=y, Z)
      fit1<-glm(trait~.,family="binomial",data=tdat1)
      mu <-fitted.values(fit1)
      y.new = y.perm - mu
      Q.perm = t(y.new) %*% K %*% y.new/2
      G.perm[i] = as.numeric(Q.perm)
    }
    perm.pval = sum(G.perm > skat.stat)/perm
  }
  else perm = "NULL"
  res = list(skat.stat = skat.stat, asym.pval = asym.pval, 
             perm.pval = perm.pval)
  return(res)
}

#bootstrap method
# pboot = function(y,G,Z,B=200){
#   t0 = TwoStep(y,G,Z)
#   t = matrix(0,B,2) 
#   tdat1<-data.frame(trait=y, Z)
#   fit1<-glm(trait~.,family="binomial",data=tdat1)
#   mu <-fitted.values(fit1)
#   X = cbind(G,Z)
#   for( i in 1:B){
#     perm= sample(1:n,n,replace = TRUE)
#     new.G = X[perm,1:m]
#     new.Z = X[perm,-(1:m)]
#     new.y=  rbinom(n,1,mu)
#     t[i,] = TwoStep(new.y,new.G,new.Z)
#   }
#   p=sapply(1:2,function(j){mean(t[,j]>t0[j])})
#   return(p)
# }

