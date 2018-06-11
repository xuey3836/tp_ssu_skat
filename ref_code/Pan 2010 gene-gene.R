###SumSqU and SumSqUw tests for main-effects and interactions for gene 1
### AND INCLUDING the main effects of gene 2 (differing from Wang et al (2009)).
### WP, 12/17/08


### Copyright: Wei Pan, 12/17/08
###            weip@biostat.umn.edu
###            http://www.biostat.umn.edu/~weip/

#input: sdat: a data-frame with components named
#              1) trait: a vector of 0's and 1's for the disease indicator;
#              2) Ro: genotypes for gene 1
#              3) Rt: genotypes for gene 2
#       B: # of simulations used to estimate the p-value for UminP; default=1000
# output: pvalues for 
#         1) main-effects + (2-way) interactions model: SumSqU, SumSqUw, score, (simulation-based) UminP, (Bonferonni-adjusted) UminP
#         2) main-effects model: SumSqU, SumSqUw, score, (simulation-based) UminP, (Bonferonni-adjusted) UminP


SumSqUsInter<-function(sdat, B=1000){
  n<-dim(sdat)[1]
  R1<-sdat[,substring(names(sdat),1,2)=="Ro"]
  mm1<-dim(R1)[2]
  R2<-sdat[,substring(names(sdat),1,2)=="Rt"]
  mm2<-dim(R2)[2]
  
  # main-effects plus interactions with gene 1:
  X<-R1
  for(i in 1:mm1)
    for(j in 1:mm2)
      X<-cbind(X, R1[,i]*R2[,j])
  # standardize X:
  X0<-cbind(X, R2)
  Xc<-apply(X, 2, function(x)(x-mean(x)) )
  X0c<-apply(X0, 2, function(x)(x-mean(x)) )
  
  Y<-sdat$trait
  #fit with gene 2's genotype:
  tdat<-data.frame(trait=sdat$trait,R2)
  tfit<-glm(trait~.,family="binomial",data=tdat)
  mu<-fitted.values(tfit)
  
  # score stat:
  U<-t(X) %*% (Y-mu) 
  
  ######Cov of the WHOLE score stat vector:
  CovS0<-t(X0c) %*% diag(mu*(1-mu)) %*% X0c
  ##cov of the score stats:
  CovS<- solve(solve(CovS0)[1:(mm1+mm1*mm2), 1:(mm1+mm1*mm2)])
  
  #unweighted one:
  Tg1<- t(U) %*% U
  ###cov of the score stats: No, not the below, which ignored the nuisance para's!
  #CovS<- t(Xc) %*% diag(mu*(1-mu)) %*% Xc
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
  
  ## score test for ONLY main-effects:
  U1<-U[1:(mm1)]
  CovS1<-CovS[1:(mm1), 1:(mm1)]
  ##gInv of CovS1:
  CovS1.edecomp<-eigen(CovS1)
  CovS1.rank<-sum(abs(CovS1.edecomp$values)> 1e-8)
  inveigen<-ifelse(abs(CovS1.edecomp$values) >1e-8, 1/CovS1.edecomp$values, 0)
  P<-solve(CovS1.edecomp$vectors)
  gInv.CovS1<-t(P) %*% diag(inveigen) %*% P
  # CovS1 %*% gInv.CovS1 %*% CovS1 == CovS1
  Tscore1<- t(U1) %*% gInv.CovS1  %*% U1
  pTscore1<-as.numeric( 1-pchisq(Tscore1, CovS1.rank) )
  
  ## SumSqUs for ONLY main-effects:
  #unweighted one:
  Tg1<- t(U1) %*% U1
  ##distr of Tg1 is sum of cr Chisq_1:
  cr<-eigen(CovS1, only.values=TRUE)$values
  ##approximate the distri by alpha Chisq_d + beta:
  alpha1<-sum(cr*cr*cr)/sum(cr*cr)
  beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
  d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
  alpha1<-as.real(alpha1)
  beta1<-as.real(beta1)
  d1<-as.real(d1)
  pSumSqU1<-(1-pchisq((Tg1-beta1)/alpha1, d1))
  #paraTg1<-c(alpha, beta, d)
  
  #weighted one:
  Tg2<- t(U1) %*%  diag(1/diag(CovS1)) %*% U1
  ##distr of Tg1 is sum of cr Chisq_1:
  cr<-eigen(CovS1 %*% diag(1/diag(CovS1)), only.values=TRUE)$values
  ##approximate the distri by alpha Chisq_d + beta:
  alpha2<-sum(cr*cr*cr)/sum(cr*cr)
  beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
  d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
  alpha2<-as.real(alpha2)
  beta2<-as.real(beta2)
  d2<-as.real(d2)
  pSumSqUw1<-(1-pchisq((Tg2-beta2)/alpha2, d2))
  #paraTg2<-c(alpha, beta, d)
  
  #### V1 %*% t(V1)=CovS:
  #CovS.edecomp<-eigen(CovS)
  CovS.ev<-ifelse(abs(CovS.edecomp$values) >1e-8, CovS.edecomp$values, 0)
  V1<-CovS.edecomp$vectors %*% diag(sqrt(CovS.ev))
  
  
  ##univariate/marginal tests: Main-effect or All:
  Tus<-as.vector(U^2/diag(CovS))
  pTus<- 1-pchisq(Tus, 1)
  pUmin<-min(pTus)
  pUminM<-min(pTus[1:(mm1)])
  
  ############################MC simulations to estimate the p-values:
  NMCsim<-B
  pUmins<-pUminMs<-rep(1, NMCsim)
  for(MCsim in 1:NMCsim){
    #########SumSqUs:
    Us<-rnorm(nrow(U))
    Us<-V1 %*% Us
    Tus<-as.vector(Us^2/diag(CovS))
    pTuss<- 1-pchisq(Tus, 1)
    pUmins[MCsim]<-min(pTuss)
    pUminMs[MCsim]<-min(pTuss[1:(mm1)])
  }
  
  return(c(pSumSqU, pSumSqUw, pTscore, 
           pUmin*length(U), sum(pUmins<pUmin)/NMCsim,
           pTscore1, pSumSqU1, pSumSqUw1,
           pUminM*(mm1), sum(pUminMs<pUminM)/NMCsim))
  
}
