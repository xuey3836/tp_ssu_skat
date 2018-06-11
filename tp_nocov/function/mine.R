
#SSU statistics
stat_ssu_skat<-function(y,X,n ,m){
  # score stat:
  y.new = y - mean(y)
  U<-t(X) %*% (y.new)
  #centering X:
  Xc<-apply(X, 2, function(x)(x-mean(x)) )
  ##cov of the score stats:
  CovS<- mean(y)*(1-mean(y))*(t(Xc) %*% Xc)
  #weighted one:
  Tg2<- t(U) %*%  diag(1/diag(CovS)) %*% U
  ##skat.stat
  # K1= kernel_IBS(2*X, n, m)
  K= IBS_kernel(X,n = n,m = m)
  
  Q = t(y.new) %*% K %*% y.new/2
  skat.stat = as.numeric(Q)
  return(list(SSU= as.numeric(Tg2),SKAT =as.numeric(skat.stat)))
}

# 
# #HWET to select the genetic model
# def.genmodel<-function(y,x){
#   tab=table(y,x)
#   r= sum(y==1)
#   s=sum(y==0)
#   n=r+s
#   colnames<-colnames(tab)
#   r1= ifelse("1" %in% colnames,tab["1","1"],0)
#   r2=ifelse("2" %in% colnames,tab["1","2"],0)
#   p1= r1/r
#   p2= r2/r
#   n1= ifelse("1" %in% colnames,colSums(tab)[2],0)
#   n2= ifelse("2" %in% colnames,colSums(tab)[3],0)
#   delta_p=p2-(p2+p1/2)^2
#   t2 = (1-n2/n-n1/(2*n))*(n2/n+n1/(2*n))
#   Z= sqrt(r)*(delta_p)/t2
#   cutpoint=qnorm(0.95)
#   if (Z>cutpoint){
#     return("R")
#   }else if(Z<(-cutpoint)){
#     return("D")
#   }else{return("A")}
# }

TwoStep<-function(y,X){
  STAT=NULL
  m=dim(X)[2]
  ##select gentic model
  gen.model=apply(X,2,function(x){
    # def.genmodel(y=y,x=x)
    defmodel(y=y,x=x)
  })
  GX=sapply(1:m,function(j){if(gen.model[j]=="R"){
    as.numeric(X[,j]>1)
  }else if(gen.model[j]=="D"){
    as.numeric(X[,j]>0)
  }else{X[,j]/2}})
  
  ## compute statistics
  stat=stat_ssu_skat(y=y,X=GX,n=length(y),m=ncol(X))
  STAT["tpSSU"] =stat$SSU
  STAT["tpSKAT"] = stat$SKAT
  return(STAT)
}

# stat.boot=boot(data=data,statistic = stat,R=2000)
# p = sapply(1:4, function(j){mean(stat.boot$t[,j]>stat.boot$t0[j])})
# B=1000,permutation
ptwostep<-function(y,X,B){
  X = as.matrix(X)
  Stat=TwoStep(y = y,X=X)
  nc=length(Stat)
  T2=matrix(0,ncol=nc,nrow = B)
  for(i in 1:B){
    perm.sample = sample(1:length(y))
    y.perm = y[perm.sample]
    T2[i,]=TwoStep(y = y.perm,X=X)
  }
  p = sapply(1:nc,function(j){mean(T2[,j]>Stat[j])})
  names(p)<-names(Stat)
  return(p)
}

## change from the AssoctesteR package:: SKAT
SKAT.m = function (y, X, perm =0) 
{
  X = as.matrix(X)
  n = nrow(X)
  m = ncol(X)
  K = IBS_kernel(X/2 ,n=n,m=m)
  y.new = y - mean(y)
  Q = t(y.new) %*% K %*% y.new/2
  mu = rep(mean(y), n)
  V = diag(mu * (1 - mu))
  ones = rep(1, n)
  P = V - V %*% ones %*% solve(t(ones) %*% V %*% ones) %*% 
    t(ones) %*% V
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
    x.perm = rep(0, perm)
    for (i in 1:perm) {
      perm.sample = sample(1:length(y))
      y.perm = y[perm.sample]
      y.new = y.perm - mean(y.perm)
      Q.perm = t(y.new) %*% K %*% y.new/2
      x.perm[i] = as.numeric(Q.perm)
    }
    perm.pval = sum(x.perm > skat.stat)/perm
  }
  else perm = "NULL"
  name = "SKAT: Sequence Kernel Association Test"
  arg.spec = c(sum(y), length(y) - sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")
  res = list(skat.stat = skat.stat, asym.pval = asym.pval, 
             perm.pval = perm.pval, args = arg.spec, name = name)
  class(res) = "assoctest"
  return(res)
}


stat4=function(y,X){
  # newdata = data[ind,]
  X = as.matrix(X)
  tp.stat=TwoStep(y,X)
  s.stat = stat_ssu_skat(y,X/2,n=nrow(X),m=ncol(X))
  c(tp.stat,s.stat$SSU,s.stat$SKAT)
}

pboot = function(y,X,B){
  t0 = stat4(y,X)
  t = matrix(0,Bt,4) 
  y_p = mean(y==1)
  for( i in 1:B){
    new.X = sample_n(as.data.frame(X),n,replace = TRUE)
    new.y=  rbinom(n,1,y_p)
    t[i,] = stat4(new.y,new.X)
  }
  p=sapply(1:4,function(j){mean(t[,j]>t0[j])})
  return(p)
}
