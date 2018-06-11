


library(mvtnorm)
library(Rcpp)
##  generate multiple correlated SNPs based on the multivariate normal distribution
## MAF.vec: the MAF vector for m SNPs; sig.mat: the covariance matrix for the multinormal distribution
## n.case and n.control are the sample size of the case and control groups; 
## bet:  the predefined coefficient used to generate case data; N: the predefined total sample size for the case data

gen.MultiSNP <- function(MAF.vec=rep(0.15,100), sig.mat, n.case, n.control,
                         beta0=0, bet=c(rep(0.2,50),rep(-0.1,50)),gamma=0.5, gen.model=rep("A",100),N=10000)
{
     m <- length(MAF.vec)
     ## calculate the cutpoint
     # the frequenticy matrix of the all SNPs (every column represent a SNP)
     fre.mat <- rbind((1-MAF.vec)^2, 2*MAF.vec*(1-MAF.vec), MAF.vec^2) 
     cum.fre <- rbind((1-MAF.vec)^2, 1-MAF.vec^2)
     cutPoint <- qnorm(cum.fre)
     
     ### Firstly, generate the all data
     dat <- rmvnorm(N, mean=rep(0,m), sigma=sig.mat) 
     SNP = sapply(1:m,  function(j){cut(x = dat[,j],cutPoint[,j])})

     ###Reccsive model
     GSNP=sapply(1:m,function(j){if(gen.model[j]=="R"){
                            as.numeric(SNP[,j]>1)
                          }else if(gen.model[j]=="D"){
                            as.numeric(SNP[,j]>0)
                          }else{SNP[,j]/2}})
     
     
     z = rnorm(N)
     effect <- exp(-(matrix(GSNP, nrow=N, ncol=m) %*% matrix(bet, ncol=1)+beta0+z*gamma))       
     y.prob <- 1/(1+effect)
     res01 <- rep(NA, N)
     for(i in 1:N)
     {
          res01[i] <- rbinom(1, 1, y.prob[i])
     }
     loc <- which(res01==1)
     loc.control<- which(res01==0)
     SNPx = cbind(SNP,z)
     SNP.case <- SNPx[sample(loc,n.case),]
     SNP.control<-SNPx[sample(loc.control,n.control),]
     data = data.frame(cbind(y=rep(c(0,1),c(n.control,n.case)),rbind(SNP.control,SNP.case)))
     
     list(SNP.control=SNP.control,SNP.case= SNP.case,df=data) 
}


cppFunction('NumericVector cut(NumericVector x, NumericVector cutpoint){
            int n = x.size();
            NumericVector g(n);
            for (int i=0;i<n;i++){
            if (x(i)<cutpoint(0)){
            g(i)=0;
            }else if(x(i)<cutpoint(1)){
            g(i)=1;
            }else{
            g(i)=2;
            }
            }
            return g;
            }')

