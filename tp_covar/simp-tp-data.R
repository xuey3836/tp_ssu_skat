# setwd("C:/Users/o0/Desktop/RR/")
rm(list = ls())
setwd("/home/o0/Desktop/tpSSU/tp_covar/")
source("pvalue-tp.R")

path="sim_data/m10_gen0R0D1A_smCS_beta0.Rdata"
index = regexpr("beta[0-9]",path)
savename= substr(path,index[1],index[1]+4)
load(path)
Rep = length(sim_data)
m = as.numeric(substr(path,11,12))
ptm=proc.time()
result=matrix(0, nrow = 4, ncol = Rep)

for (r in 1:Rep){
  if (r%%100==0){
    cat("the loop is", r ,"\n")
  }
  ##generate data
  p=NULL
  dat = sim_data[[r]]
  y=dat$y
  G = as.matrix(dat[,2:(m+1)])
  Z = as.matrix(dat$z)
  result[,r]=pvalue_tp(y,G,Z,B=200)
}

proc.time()-ptm
power = apply(result<0.05,1,mean)
path
power
write.csv(t(power),paste0("result/",savename,".csv"))



