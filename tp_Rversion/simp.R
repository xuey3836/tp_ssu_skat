 setwd("C:/Users/liuyan/Desktop/case-only/")
rm(list = ls())
#setwd("/home/o0/Desktop/case-only/")
source("function/gen_MultiSNP.r")
source("function/mine.R")
###load package
library(mvtnorm)
library(AssotesteR)

#set parameters
m=10

####Sig.mat
sig<-c("CS","AR-1","Rand")
sig.index=sig[1]
##AR-1
sig.mat=matrix(0,m,m)
for (i in 1:m){
  for(j in i:m){
    if(sig.index=="CS"){
      sig.mat[i,j]=0.4
    }else if(sig.index=="AR-1"){
      sig.mat[i,j]=0.8^abs(i-j)
    }else if(sig.index=="Rand"){
      sig.mat[i,j]=runif(1,0.3,0.7)
    }
    sig.mat[j,i] = sig.mat[i,j]
  }
}
diag(sig.mat)=1

#disease-casusing SNP
nlociR=0
nlociD=0
lociR=sample(1:m,nlociR,replace = FALSE)
lociD=sample(setdiff(1:m,lociR),nlociD,replace = FALSE)

gen.model=rep("A",m)
gen.model[lociR]=rep("R",nlociR)
gen.model[lociD]=rep("D",nlociD)

bet=rep(0,m)
beta0=log(0.05/0.95)
# beta.v=log(c(1,1.2,1.4,1.6,1.8,2))
beta.v=log(seq(1,3,0.2))
lociA=NULL
nlociA=0
if(nlociD+nlociR==0){
  nlociA=1
  beA.index= 1
  lociA=sample(1:m,nlociA,replace = FALSE)
  bet[lociA]=beta.v[beA.index]
}

nloci=nlociA+nlociR+nlociD
loci=c(lociA,lociR,lociD)
##MaF
MAF.vec=runif(m,0.2,0.8)
MAF.vec[loci]=rep(0.3,nloci)
#Beta

beR.index =1
beD.index= 1
##null 
bet[lociR[1]]=beta.v[beR.index]
bet[lociR[2]]=-beta.v[beR.index]
bet[lociD]=beta.v[beD.index]

n.case=500
n.control=500
n=n.case+n.control
N=300000

Rep=1
Bt=1000
set.seed(12345)
seed = sample(1:2000000,Rep,replace = FALSE)
ptm=proc.time()
result=matrix(0, nrow = 6, ncol = Rep)

for (r in 1:Rep){
  set.seed(seed[r])
  if (r%%100==0){
    cat("the loop is", r ,"\n")
  }
  ##generate data
  p=NULL
  set.seed(11)
  SNP =gen.MultiSNP(MAF.vec = MAF.vec,sig.mat = sig.mat,n.case=n.case,
                    beta0 = beta0,n.control = n.control,N=N,bet = bet,
                    gen.model=gen.model)
  data = SNP$df
  y=data$y
  X=as.matrix(data[,-1])
  ##two-SSU
  p = ptwostep(y=y,X=X,B=Bt)
  
  # SSU and SKAT
  ssuw = SSUW(y=y,X = X,perm=Bt)
  p["SSUa"]=ssuw$asym.pval
  p["SSUp"]= ssuw$perm.pval
  skat = SKAT(y=y,X=X,kernel = "IBS",perm = Bt)
  p["SKATa"]=skat$asym.pval
  p["SKATp"] = skat$perm.pval
  result[,r]=p
}

proc.time()-ptm
power=apply(result<0.05,1,mean)
power
write.csv(t(power),paste0("result/","power",beR.index,".csv"))

# 
# ##arrange the result
# power = NULL
# for (l in 2:6){
#   p = read.csv(file=paste0("result/","power",l,".csv"))
#   p[1] = l
#   power = rbind(power,p)
# }
# write.csv(power,"power.csv")





