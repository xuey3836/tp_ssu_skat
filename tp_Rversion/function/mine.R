
#SSU statistics
SSU.stat<-function(y,X,gen.model){
  # score stat:
  U<-t(X) %*% (y-mean(y))
  
  #centering X:
  Xc<-apply(X, 2, function(x)(x-mean(x)) )
  #unweighted one:
 # Tg1<- t(U) %*% U
  ##cov of the score stats:
  CovS<- mean(y)*(1-mean(y))*(t(Xc) %*% Xc)
  #weighted one:
  Tg2<- t(U) %*%  diag(1/diag(CovS)) %*% U
  #return(list(SSU.stat=as.numeric(Tg1),SSUw.stat=as.numeric(Tg2)))
  return(list(SSUw.stat=as.numeric(Tg2)))
  
  
}

##SKAT statistics
SKAT.stat=function (y, X, kernel = "linear", weights = NULL, a = 0.5, b = 0.5, 
          perm = NULL) 
{
  n = nrow(X)
  p = ncol(X)
 
  K = switch(kernel, linear = X %*% t(X), 
             wlinear = ((X %*%  diag(weights^2)) %*% t(X)), 
             quadratic = (X %*% t(X) +  1)^2, 
             IBS = kernel_IBS(X, n, p))
  y.new = y - mean(y)
  Q = t(y.new) %*% K %*% y.new/2
  skat.stat = as.numeric(Q)
  return(skat.stat)
}

kernel_IBS <-
  function(Z, n, p)
  {
    K = diag(1, n, n)	
    aux = .C("kernel_IBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)))[[4]]
    matrix(aux, nrow=n)
  }

##HWET to select the genetic model
def.genmodel<-function(y,x){
  tab=table(y,x)
  r= sum(y==1)
  s=sum(y==0)
  n=r+s
  colnames<-colnames(tab)
  r1= ifelse("1" %in% colnames,tab["1","1"],0)
  r2=ifelse("2" %in% colnames,tab["1","2"],0)
  p1= r1/r
  p2= r2/r
  n1= ifelse("1" %in% colnames,colSums(tab)[2],0)
  n2= ifelse("2" %in% colnames,colSums(tab)[3],0)
  delta_p=p2-(p2+p1/2)^2
  t2 = (1-n2/n-n1/(2*n))*(n2/n+n1/(2*n))
  Z= sqrt(r)*(delta_p)/t2
  cutpoint=qnorm(0.95)
  if (Z>cutpoint){
    return("R")
  }else if(Z<(-cutpoint)){
    return("D")
  }else{return("A")}
}


TwoStep<-function(y,X){
  STAT=NULL
  m=dim(X)[2]
  ##select gentic model
  gen.model=apply(X,2,function(x){
    def.genmodel(y=y,x=x)
  })
  GX=sapply(1:m,function(j){if(gen.model[j]=="R"){
    as.numeric(X[,j]>1)
  }else if(gen.model[j]=="D"){
    as.numeric(X[,j]>0)
  }else{X[,j]/2}})
  
  ## compute statistics
  ssu=SSU.stat(y=y,X=GX)
  STAT["tpSSU"]=ssu$SSUw.stat
  GX2 = 2*GX
  STAT["tpSKAT"]=SKAT.stat(y=y,X=GX2,kernel = "IBS")
  
  return(STAT)
}

# B=1000
ptwostep<-function(y,X,B){
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






