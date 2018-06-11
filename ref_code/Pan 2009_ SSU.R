###SumSqU, SumSqUw, UminP, score  tests  for input (Y, X);
### X may contain 1) main-effects and/or
### 2) main-ffects+interactions.

### Copyright: Wei Pan, 1/3/09
###            weip@biostat.umn.edu 
###            http://www.biostat.umn.edu/~weip/

#input: Y: a vector of 0's and 1's for the disease indicator
#       X: a genotype matrix; each row is for one subject, and each col for 
#          one SNP (with 1-DF coding by additive, recessive or dominant 
#          inheritance mode) or one of the two indicators for the 2-DF coding; 
#          for additive mode, coded as 0, 1 and 2.
#       B: # of simulations used to estimate the p-value for UminP; default=1000
# output: pvalues for SumSqU, SumSqUw, minP without adjustment for multiple 
#         testing (hence NOT valid), UminP and score tests, respectively.

SumSqUs<-function(Y, X, B=1000){
  # score stat:
  U<-t(X) %*% (Y-mean(Y))
  
  #centering X:
  Xc<-apply(X, 2, function(x)(x-mean(x)) )
  
  #unweighted one:
  Tg1<- t(U) %*% U
  ##cov of the score stats:
  CovS<- mean(Y)*(1-mean(Y))*(t(Xc) %*% Xc)
  ##distr of Tg1 is sum of cr Chisq_1:
  cr<-eigen(CovS, only.values=TRUE)$values
  ##approximate the distri by alpha Chisq_d + beta:
  alpha1<-sum(cr*cr*cr)/sum(cr*cr)
  beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
  d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
  alpha1<-as.real(alpha1)
  beta1<-as.real(beta1)
  d1<-as.real(d1)
  pSumSqU<-(1-pchisq((Tg1-beta1)/alpha1, d1))
  #paraTg1<-c(alpha, beta, d)
  
  #weighted one:
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
  pSumSqUw<-(1-pchisq((Tg2-beta2)/alpha2, d2))
  #paraTg2<-c(alpha, beta, d)
  
  ##standard score test:
  ##gInv of CovS:
  CovS.edecomp<-eigen(CovS)
  CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
  inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
  P<-solve(CovS.edecomp$vectors)
  gInv.CovS<-t(P) %*% diag(inveigen) %*% P
  # CovS %*% gInv.CovS %*% CovS == CovS
  Tscore<- t(U) %*% gInv.CovS  %*% U
  pTscore<-as.numeric( 1-pchisq(Tscore, CovS.rank) )
  
  #### V1 %*% t(V1)=CovS:
  #CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(CovS.ev))
  
  ##univariate/marginal tests: 
  Tus<-as.vector(U^2/diag(CovS))
  pTus<- 1-pchisq(Tus, 1)
  pUmin<-min(pTus)
  ############################MC simulations to estimate the p-values:
  NMCsim<-B
  pUmins<-rep(1, NMCsim)
  for(MCsim in 1:NMCsim){
    #########SumSqUs:
    Us<-rnorm(nrow(U))
    Us<-V1 %*% Us
    Tus<-as.vector(Us^2/diag(CovS))
    pTuss<- 1-pchisq(Tus, 1)
    pUmins[MCsim]<-min(pTuss)
  }
  pTu<-sum(pUmins<pUmin)/NMCsim
  pTuB<-pUmin*length(U)
  
  c(pSumSqU, pSumSqUw, pTu, pTuB, pTscore)
}