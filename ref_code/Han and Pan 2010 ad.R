######################Adaptive Sum Test
######################modified from the program of Fang HAN (Email: fanghan@biostat.umn.edu) on 1/22/10

library("mvtnorm")

# Introduction to input parameters:

# Y: disease lables; =0 for controls, =1 for cases;
# X: genotype data; each row for a subject, and each column for an SNP;
#    the value of each element is the # of the copies for an allele;
# alpha0 = 0.1: the cut-off value; if the marginal p-value of an SNP is less 
#               than alpha0 and its marginal regression coefficient is negative,
#               we flip its coding;
# B=1000: the number of permutations; 

############################
# Output is the p-values of the seven tests in the order of:                 
# 1. Score 
# 2. SSU
# 3. SSUw
# 4. UminP
# 5. Sum
# 6. aSum-P
# 7. aSum


############################################
#######Some back-up programs
############################################

PowerUniv <- function(U,V){
  n <- dim(V)[1]
  
  x <- as.numeric(max(abs(U)))
  TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
  
  return(TER)
}

SumTest <- function(Y,X,alpha0){
  pv<-NULL
  beta<-NULL
  for(i in 1:ncol(X)){
    fit  <- glm(Y~X[,i], family=binomial(logit))
    beta <- c(beta,fit$coefficients[-1])
    pv <- c(pv,as.numeric(summary(fit)$coefficients[,4][-1]))
  }
  
  Xg <- X
  
  Xbar<-mean(Xg) 
  Xgb<-Xg
  
  Xgb<-Xgb-Xbar
  
  U<-t(Xg) %*% (Y-mean(Y))
  CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
  
  a<-rep(1, length(U))
  a[beta<0 & pv<alpha0] <- -1
  
  u <- sum(t(a)%*%U)
  v <- as.numeric(t(a) %*% CovS %*% (a))
  
  Tsum<- sum(t(a)%*%U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
  pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )       
  
  return(cbind(u,pTsum,v,a*U,a))
}

rSumTest <- function(Y,X,k,alpha0){
  nCase <- sum(Y==1)
  nControl <- sum(Y==0)
  n <- nCase+nControl
  
  u0<-pv<-v0<-u00<-NULL    
  for(sim in 1:k){
    set.seed(sim)
    case <- sample(1:n,nCase,replace = FALSE)
    Y[case] <- 1
    Y[-case] <- 0   
    fit <- SumTest(Y,X,alpha0)
    u0 <- c(u0,fit[1,1])
    pv <- c(pv,as.numeric(fit[1,2]))
    v0 <- c(v0,fit[1,3])
    u00 <- cbind(u00,fit[,4])
  }   
  
  a <- rep(1,nrow(u00))  
  
  mu <- mean(u0)
  v <- var(u0)
  
  x<-(u0-mu)^2/v
  a <- sqrt(var(x)/2)
  b <- mean(x)-a
  
  return(cbind(rep(mean(u0),k),rep(v,k),rep(a,k),rep(b,k),pv))
}
############################################
############################################
############################################



aSumTest<-function(Y, X, B=1000, alpha0=0.1){
  
  Xg <- X
  
  ############SSU#####################
  Xbar<-apply(Xg, 2, mean) 
  Xgb<-Xg
  for(i in 1:nrow(Xg))
    Xgb[i,]<-Xg[i,]-Xbar
  
  #########SumSqUs:
  U<-t(Xg) %*% (Y-mean(Y))
  
  #SumSqU:
  Tg1<- t(U) %*% U
  ##cov of the score stats:
  CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
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
  nSNP<-nrow(CovS)  # added by WP, 10/24/10
  Vs <- matrix(c(rep(0,nSNP^2)),nrow=nSNP)
  for(i in 1:nSNP){
    for(j in 1:nSNP){
      Vs[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
    }
  }
  pTus <- as.numeric(PowerUniv(Tus,Vs))
  
  
  
  ##########SumTest########################
  a<-rep(1, length(U))
  Tsum<- sum(t(a)%*%U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
  pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )  
  
  
  ##########aSumTest
  fit0 <- SumTest(Y,Xg,alpha0)
  u <- fit0[1,1]
  pv <- fit0[1,2]
  v <- fit0[1,3]
  fit <- rSumTest(Y,Xg,B,alpha0)
  u0 <- fit[1,1] 
  v0 <- fit[1,2]
  a <- fit[1,3]
  b <- fit[1,4]
  pv0 <- fit[,5]
  
  aSumP <- sum(pv>pv0)/length(pv0)
  aSum <- as.numeric( 1-pchisq(abs(((u-u0)^2/v0-b)/a), 1) )
  
  return(c(pTscore, pTg1, pTg2, pTus, pTsum, aSumP, aSum))
}  


