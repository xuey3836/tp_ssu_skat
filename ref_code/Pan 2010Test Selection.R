# April 6, 2010 
# Programmer: Fang HAN
# Email: fanghan@biostat.umn.edu



####################################
##############Base Function ########
####################################


#################Function of Computing different test power ########################
library(mvtnorm)
#alpha <- 0.05


TrunPlus <- function(x,c){
  y <- x-c
  return(y*as.numeric(y>0)+c)
}


########################### The Score Statistics ###############################
PowerScore <- function(U,V,alpha){
  n <- dim(V)[1]
  I <- matrix(c(rep(0,n^2)),nrow=n)
  I[row(I)==col(I)]<-1
  c <- qchisq(1-alpha, df=n)
  
  CovS<-V
  ##gInv of CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
  inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
  P<-solve(CovS.edecomp$vectors)
  gInv.CovS<-t(P) %*% diag(inveigen) %*% P
  
  TER <- 1 - pchisq(c, df=n, ncp = (t(U)%*%gInv.CovS%*%U)[1])
  return(TER)
}

PowerScoreBias <- function(V,N,alpha){
  n <- dim(V)[1]
  
  mu0 <- array(c(rep(0,n)))
  Sum <- 0
  
  CovS <- V
  #### V1 %*% t(V1)=CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(abs(CovS.ev)))
  
  
  for(i in 1:N){
    Us<-rnorm(n)
    Us<-V1 %*% Us
    
    Sum <- Sum + PowerScore(array(Us),V,alpha)
  }
  return(TrunPlus(Sum/N-alpha,0))
}

PowerScoreEst <- function(U,V,alpha){
  return(PowerScore(array(U),V,alpha))
}



########################### Univariate Statistics ###############################
PowerUniv <- function(U,V,alpha){
  n <- dim(V)[1]
  
  x <- as.numeric(max(abs(U)))
  TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
  
  return(TER)
}


rPowerUniv <- function(Tus,Vs,c,alpha){
  n <- length(Tus)        
  cdfreal <- function(x){
    pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=as.numeric(array(Tus)),sigma=Vs)         
  }
  
  return(1 - cdfreal(c))
}



PowerUnivBias <- function(Vs,c,N,alpha){
  n <- dim(Vs)[1]
  
  mu0 <- array(c(rep(0,n)))
  Sum <- 0
  
  CovS <- Vs
  #### V1 %*% t(V1)=CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(abs(CovS.ev)))
  
  for(i in 1:N){
    Us<-rnorm(n)
    Us<-V1 %*% Us
    
    cdfreal <- function(x){
      pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=as.numeric(array(Us)),sigma=Vs)         
    }
    
    Sum <- Sum + 1 - cdfreal(c)
  }
  
  return(TrunPlus(Sum/N-alpha,0))
}



########################### SumSqU Statistics #############################
PowerSumSqU <- function(U,V,alpha){
  n <- dim(V)[1]
  c1 <- eigen(V)$values
  c1sum <- sum(c1)
  c1s <- sum(c1^2)
  c1t <- sum(c1^3)
  
  a1 <- c1t/c1s
  b1 <- c1sum-c1s^2/c1t
  d1 <- c1s^3/c1t^2
  
  a1 <- as.real(a1)
  b1 <- as.real(b1)
  d1 <- as.real(d1)
  
  c <- qchisq(1-alpha, df=d1)
  
  TER <- 1 - pchisq(c, df=d1, ncp = sum(U^2)/a1)
  
  return(TER)    											
}

PowerSumSqUBias <- function(V,N,alpha){
  n <- dim(V)[1]
  
  mu0 <- array(c(rep(0,n)))
  Sum <- 0
  
  CovS <- V
  #### V1 %*% t(V1)=CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(abs(CovS.ev)))
  
  
  for(i in 1:N){
    Us<-rnorm(n)
    Us<-V1 %*% Us
    
    Sum <- Sum + PowerSumSqU(array(Us),V,alpha)
  }
  return(TrunPlus(Sum/N-alpha,0))
}

PowerSumSqUEst <- function(U,V,alpha){
  return(PowerSumSqU(array(U),V,alpha))
}




########################### SumSqUw Statistics #############################
PowerSumSqUw <- function(U,V,alpha){
  n <- dim(V)[1]
  
  c1 <- eigen(V%*%diag(1/diag(V)))$values
  c1sum <- sum(c1)
  c1s <- sum(c1^2)
  c1t <- sum(c1^3)
  
  a1 <- c1t/c1s
  b1 <- c1sum-c1s^2/c1t
  d1 <- c1s^3/c1t^2
  
  a1 <- as.real(a1)
  b1 <- as.real(b1)
  d1 <- as.real(d1)
  
  c <- qchisq(1-alpha, df=d1)
  
  TER <- 1 - pchisq(c, df=d1, ncp = (t(U)%*%diag(1/diag(V))%*%U)[1]/a1)
  
  return(TER)    											
}

PowerSumSqUwBias <- function(V,N,alpha){
  n <- dim(V)[1]
  
  mu0 <- array(c(rep(0,n)))
  Sum <- 0
  
  CovS <- V
  #### V1 %*% t(V1)=CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(abs(CovS.ev)))
  
  
  for(i in 1:N){
    Us<-rnorm(n)
    Us<-V1 %*% Us
    
    Sum <- Sum + PowerSumSqUw(array(Us),V,alpha)
  }
  return(TrunPlus(Sum/N-alpha,0))
}

PowerSumSqUwEst <- function(U,V,alpha){
  return(PowerSumSqUw(array(U),V,alpha))
}





########################### The Sum Statistics ###############################
PowerSum <- function(U,V,alpha){
  n <- dim(V)[1]
  One <- c(rep(1,n))
  Uc <- t(One)%*%U
  Vc <- t(One)%*%V%*%One   
  Uc2 <- Uc/sqrt(Vc)
  
  c <- qchisq(1-alpha, df=1)
  
  TER <- 1 - pchisq(c, df=1, ncp = t(Uc2)%*%Uc2)
  
  return(TER)
}

PowerSumBias <- function(V,N,alpha){
  n <- dim(V)[1]
  
  mu0 <- array(c(rep(0,n)))
  Sum <- 0
  
  CovS <- V
  #### V1 %*% t(V1)=CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(abs(CovS.ev)))
  
  
  for(i in 1:N){
    Us<-rnorm(n)
    Us<-V1 %*% Us
    
    Sum <- Sum + PowerSum(array(Us),V,alpha)
  }
  return(TrunPlus(Sum/N-alpha,0))
}

PowerSumEst <- function(U,V,alpha){
  return(PowerSum(array(U),V,alpha))
}



TheoryC <- function(U,V,rtol,alpha){
  n <- dim(V)[1]
  
  cdf <- function(x){
    pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V)
  }
  
  target<-function(x){
    cdf(x)-(1-alpha) 
  }
  
  cwhole <- uniroot(target, c(0, 3*max(V)), tol = rtol)
  c<-cwhole$root
  
  return(c)
}




PowerMax <- function(U,V,N,bias1,bias2,bias3,bias4,bias5,c,V1,M,tm,alpha){
  n <- dim(V)[1]
  nSNP<-n
  mu <- c(rep(0,n)) 
  
  CovS <- V
  
  pTS<-c(rep(0,N)) 
  pCmin<-pCFisher<-pCTPM<-c(rep(0,N))
  
  for(sim in 1:N){
    Us<-rnorm(n)
    Us<-V1 %*% Us
    
    
    #SumSqU:
    Tg1<- t(Us) %*% Us
    ##cov of the score stats:
    CovS<- V
    ##distr of Tg1 is sum of cr Chisq_1:
    cr<-eigen(CovS, only.values=TRUE)$values
    ##approximate the distri by alpha Chisq_d + beta:
    alpha1<-sum(cr*cr*cr)/sum(cr*cr)
    beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
    d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
    alpha1<-as.real(alpha1)
    beta1<-as.real(beta1)
    d1<-as.real(d1)
    pTg1<-1-pchisq((Tg1-beta1)/alpha1, d1)
    
    
    #SumSqUw:
    Tg2<- t(Us) %*%  diag(1/diag(CovS)) %*% Us
    ##distr of Tg1 is sum of cr Chisq_1:
    cr<-eigen(CovS %*% diag(1/diag(CovS)), only.values=TRUE)$values
    ##approximate the distri by alpha Chisq_d + beta:
    alpha2<-sum(cr*cr*cr)/sum(cr*cr)
    beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
    d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
    alpha2<-as.real(alpha2)
    beta2<-as.real(beta2)
    d2<-as.real(d2)
    pTg2<-1-pchisq((Tg2-beta2)/alpha2, d2)
    
    
    ##########score test:
    ##gInv of CovS:
    CovS.edecomp<-eigen(CovS)
    CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
    inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
    P<-solve(CovS.edecomp$vectors)
    gInv.CovS<-t(P) %*% diag(inveigen) %*% P
    Tscore<- t(Us) %*% gInv.CovS  %*% Us
    pTscore<-as.numeric( 1-pchisq(Tscore, CovS.rank) )
    
    #### V1 %*% t(V1)=CovS:
    CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
    V1<-CovS.edecomp$vectors %*% diag(sqrt(CovS.ev))
    
    ##univariate/marginal tests:
    
    Tus<-as.vector(abs(Us)/sqrt(diag(CovS)))
    Vs <- matrix(c(rep(0,nSNP^2)),nrow=nSNP)
    for(i in 1:nSNP){
      for(j in 1:nSNP){
        Vs[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
      }
    }
    #c <- TheoryC(Tus,Vs)
    pTus <- PowerUniv(Tus,Vs)
    
    
    ##Sum test:
    a<-rep(1, length(Us))
    Tsum<- sum(Us)/sqrt(as.numeric(t(a) %*% CovS %*% (a)))
    pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
    
    ##PowCom test:
    pp1 <- PowerScoreEst(Us,CovS,alpha)-bias1
    pp2 <- rPowerUniv(Tus,Vs,c,alpha)-bias2
    pp3 <- PowerSumSqUEst(Us,CovS,alpha)-bias3
    pp4 <- PowerSumSqUwEst(U,CovS,alpha)-bias4
    pp5 <- PowerSumEst(U,CovS,alpha)-bias5
    
    ######## Propose Four Combined P #######
    prank <- rank(c(pp1,pp2,pp3,pp4,pp5),ties.method= tm)
    
    pTS[sim] <- pTscore^(prank[1]>4) * pTus^(prank[2]>4) * pTg1^(prank[3]>4) * pTg2^(prank[4]>4) * pTsum^(prank[5]>4)
    
    pVs<- c(pTscore,pTus,pTg1,pTg2,pTsum)
    pCmin[sim]<-min(pVs)
    pCFisher[sim]<-prod(pVs)
    pCTPM[sim]<-prod(pVs[pVs<0.05])
    
  }
  return(cbind(pTS,pCmin,pCFisher,pCTPM))
}



##################################################
##################################################




##############Introduction to function parameters

####  Y: indicator of disease status
####  X: SNP codings
####  alpha: Significance level
####  M: simulation times to derive Monte Carlo bias estimate
####  Nsim: the permutation times to estimate null hypothesis baseline p-value.
####  tol: tolerance number
####  ties.method = c("first", "random")

PowerTSTC<-function(Y,X,alpha=0.05,M=10,Nsim=1000,tol=10^(-5),tm="first"){
  
  
  nSNP <- dim(X)[2]
  #######################Achieving Score Value######################
  
  Xbar<-apply(X, 2, mean) 
  Xgb<-X
  for(i in 1:nrow(X))
    Xgb[i,]<-X[i,]-Xbar
  U<-t(Xgb) %*% (Y-mean(Y))
  CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb) 
  
  ######################applying various of test####################
  
  
  #SumSqU:
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
  
  
  #SumSqUw:
  Tg2<- t(U) %*%  diag(1/diag(CovS)) %*% U
  ##distr of Tg1 is sum of cr Chisq_1:
  cr<-eigen(CovS %*% diag(1/diag(CovS)), only.values=TRUE)$values
  ##approximate the distri by alpha Chisq_d + beta:
  alpha2<-sum(cr*cr*cr)/sum(cr*cr)
  beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
  d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
  alpha2<-as.real(alpha2)
  beta2<-as.real(beta2)
  d2<-as.real(d2)
  pTg2<-as.numeric(1-pchisq((Tg2-beta2)/alpha2, d2))
  
  
  ##########score test:
  ##gInv of CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
  inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
  P<-solve(CovS.edecomp$vectors)
  gInv.CovS<-t(P) %*% diag(inveigen) %*% P
  Tscore<- t(U) %*% gInv.CovS  %*% U
  pTscore<-as.numeric( 1-pchisq(Tscore, CovS.rank) )
  
  #### V1 %*% t(V1)=CovS:
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(CovS.ev))
  
  ##univariate/marginal tests:
  
  Tus<-as.vector(abs(U)/sqrt(diag(CovS)))
  Vs <- matrix(c(rep(0,nSNP^2)),nrow=nSNP)
  for(i in 1:nSNP){
    for(j in 1:nSNP){
      Vs[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
    }
  }
  c <- TheoryC(Tus,Vs,tol,alpha)
  pTus <- as.numeric(PowerUniv(Tus,Vs,alpha))
  
  
  ##Sum test:
  a<-rep(1, length(U))
  Tsum<- sum(U)/sqrt(as.numeric(t(a) %*% CovS %*% (a)))
  pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
  
  
  ##PowCom test:
  bias1 <- PowerScoreBias(CovS,M,alpha)
  bias2 <- PowerUnivBias(Vs,c,M,alpha)
  bias3 <- PowerSumSqUBias(CovS,M,alpha)
  bias4 <- PowerSumSqUwBias(CovS,M,alpha)
  bias5 <- PowerSumBias(CovS,M,alpha)
  
  pp1 <- PowerScoreEst(U,CovS,alpha)-bias1
  pp2 <- rPowerUniv(Tus,Vs,c,alpha)-bias2
  pp3 <- PowerSumSqUEst(U,CovS,alpha)-bias3
  pp4 <- PowerSumSqUwEst(U,CovS,alpha)-bias4
  pp5 <- PowerSumEst(U,CovS,alpha)-bias5
  
  
  
  ######## Propose one test selection and three test combinition methods #######
  prank <- rank(c(pp1,pp2,pp3,pp4,pp5),ties.method= tm)
  
  pTS <- pTscore^(prank[1]>4) * pTus^(prank[2]>4) * pTg1^(prank[3]>4) * pTg2^(prank[4]>4) * pTsum^(prank[5]>4)
  
  pVs<- c(pTscore,pTus,pTg1,pTg2,pTsum)
  pCmin<-min(pVs)
  pCFisher<-prod(pVs)
  pCTPM<-prod(pVs[pVs<0.05])
  
  ######getting Null hypothesis baseline p-value
  pmax <- PowerMax(U,CovS,Nsim,bias1,bias2,bias3,bias4,bias5,c,V1,M,tm,alpha)
  
  ######deriving four adjusted p-values
  epTS<-sum(pTS>pmax[,1])/Nsim
  epCmin<- sum(pCmin>pmax[,2])/Nsim
  epCFisher<- sum(pCFisher>pmax[,3])/Nsim
  epCTPM<- sum(pCTPM>pmax[,4])/Nsim 
  
  ######output the result
  return(cbind(pTscore,pTus,pTg1,pTg2,pTsum,epTS,epCmin,epCFisher,epCTPM)) 
}

# Output is the pvalues of the methods, the order is:
# 1: Score Test 
# 2: UminP Test
# 3: SSU Test
# 4: SSUw Test
# 5: Sum Test
# 6: Test Selection Method
# 7-9: Test combinition methods: minP, Fisher and TPM 