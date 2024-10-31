#Load packages
library("magic")#install.packages("magic")
library("mixsmsn")
library(mnormt) #install.packages("mnormt")

#Load R functions


source("densSMSN.R")
source("SMSN.mmixregre.R")
source("genmixsmsn.R")
source("utils.multivariateSMSN.R")

#Configuration
replicate   <- 300
p           <- 2
family      <- "Skew.t"


if(family == "Skew.t" || family == "Skew.slash")
{
  ## Configuracion de objetos para familia Skew.t o Skew.slash
  estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3*2 + 2*2 + 2 + 1)))
  colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Beta2.3","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2","nu")
  estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
  
  EPSMSN_1              <- array(0,dim=c(replicate,ncol=(3 + 3 + 2*2 + 2)))
  colnames(EPSMSN_1)    <- c("EPBeta1.1","EPBeta1.2","EPBeta1.3","EPBeta2.1","EPBeta2.2","EPBeta2.3","EPShape1.1","EPShape1.2","EPShape2.1","EpShape2.2","EPpii1","EPnu")
  
  EPSMSN_2              <- array(0,dim=c(replicate*p,ncol=4))
  
}

if(family == "Skew.cn")
{
  ## Configuracion de objetos para familia Skew.t o Skew.slash
  estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3*2 + 2*2 + 2 + 1 + 1)))
  colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Beta2.3","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2","nu","gama")
  estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
}

if(family == "Skew.normal")
{
  ## Configuracion de objetos para familia Skew.normal (no tiene grados de libertad)
  estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3*2 + 2*2 + 2)))
  colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Beta2.3","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2")
  estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
  EP                    <- array(0,dim=c(replicate,ncol=17))
  colnames(EP)          <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Beta2.3","Shape1.1","Shape1.2","Shape2.1","Shape2.2","sigma11.1","sigma12.1","sigma22.1","sigma11.2","sigma12.2","sigma22.2", "pii1")
}

##################################################################################################################
## Gerando dados mixturas de regressoes smsn multi
##################################################################################################################

pii    <- c(0.3,0.7)
n      <- 2000#Sample size
x      <- cbind(rep(1,n*p),runif(n*p,0,1),runif(n*p, 2,4))
#x  <- head(kronecker(cbind(rep(1,n),ais$Ferr/100),matrix(1,2,1)))
#x  <- kronecker(cbind(rep(1,n),runif(n,0,1),runif(n, 2,4)),matrix(1,2,1))


#Freedom degrees
if(family == "Skew.t" || family == "Skew.slash")
{ nu   <- 3   #No rodar para Skew.normal e normal
}else{
  nu   <- c(0.1,0.1) #Para Skew.cn
}


if(family == "Skew.t")      {k1  <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2);b <- -sqrt(2/pi)*k1}
if(family == "Skew.normal") {k1  <- 1;b <- -sqrt(2/pi)*k1}
if(family == "Skew.cn")     {k1  <- nu[1]/nu[2]^(1/2)+1-nu[1];b <- -sqrt(2/pi)*k1}
if(family == "Skew.slash")  {k1  <- 2*nu/(2*nu-1);b <- -sqrt(2/pi)*k1}

#First component
betas1 <- c(-1,-4,-3) #betas1=matrix(c(2,-3,1,7), 2,2)
Sigma1 <- matrix(c(3,1,1,3), 2,2)
shape1 <- matrix(c(1,2),2,1)
delta1 <- shape1/sqrt(1 + sum(shape1^2))
Delta1 <- matrix.sqrt(Sigma1)%*%delta1
mu1    <- t(matrix(x%*%betas1 + rep(b*Delta1,n),p,n))

#Second component
betas2 <- c(3,7,2)  #betas2=matrix(c(-1,4,-2,-4), 2,2)
Sigma2 <- matrix(c(2,1,1,2), 2,2)
shape2 <- matrix(c(1,1),2,1)
delta2 <- shape2/(sqrt(1 + sum(shape2^2)))
Delta2 <- matrix.sqrt(Sigma2)%*%delta2
mu2    <- t(matrix(x%*%betas2+rep(b*Delta2,n),p,n))

#beta      <- list(betas1,betas2)
#beta[[1]] <- matrix(c(2,-3,1,7), 2,2); beta[[2]] <- matrix(c(-1,4,-2,-4), 2,2) # matrix(c(-2,-2,-2,-3), 2,2)
betas   <- cbind(betas1,betas2)
########################################################################################################################################
j       <- 1
start.timeI      <- proc.time(); #Begin time
while(j <= replicate)
{
  cat("#Replicate",j,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')
  
  if(family == "Skew.t" || family == "Skew.slash" || family == "Skew.cn")
  {     y <- genmixsmsn(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu);hist(y[,1])
  }else{y <- genmixsmsn(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu=NULL);hist(y)}
  
  initial.values  <- initial.Values.fm.smsn.mr(y,x,g=2,get.init="k-means",family=family,lower=1.1,upper=20,space=1,plotLog = TRUE,searchNU=TRUE,printNU=FALSE, saveFigure = FALSE)
  TrySMSN         <- try(smsnReg.mmix(y, x, nu=initial.values$nu,  betas = initial.values$beta, Sigma = initial.values$Sigma, shape = initial.values$shape, pii = initial.values$pii, g = 2, get.init = FALSE, criteria = TRUE,group = FALSE, family = family, error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,obs.prob= FALSE, kmeans.param = NULL))
  #TrySMSN         <- try(smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,group = FALSE, family = "Skew.normal", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,obs.prob= FALSE, kmeans.param = NULL))
  EP              <- imm.smsn.mmrm(y, x, TrySMSN)$EP
  
  
  if(family == "Skew.t" || family == "Skew.slash" || family == "Skew.cn")
  if(class(TrySMSN)!="try-error")
    {
      print(TrySMSN$Sigma)
      estimSMSN_1[j,]             <- as.matrix(c(TrySMSN$beta[,1],TrySMSN$beta[,2],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$pii,TrySMSN$nu),byrow=TRUE)
      estimSMSN_2[(2*j-1):(2*j),] <- cbind(matrix(unlist(TrySMSN$Sigma[1]), ncol = 2, byrow = TRUE),matrix(unlist(TrySMSN$Sigma[2]), ncol = 2, byrow = TRUE))
      
      EP                          <- imm.smsn.mmrm(y, x, TrySMSN)$EP
      EPSMSN_1[j,]                <- round(c(EP[1:3],EP[9:11],EP[4:5],EP[12:13],EP[17],EP[18]),digits=4)
      print(EPSMSN_1[j,])
      EPSMSN_2[(2*j-1):(2*j),]    <- cbind(matrix(c(EP[6],EP[7],EP[7],EP[8]), ncol = 2, byrow = TRUE),matrix(c(EP[13],EP[14],EP[14],EP[15]), ncol = 2, byrow = TRUE))
      if(j>1) print(rbind(round(apply(EPSMSN_1[1:j,],2,mean),digits = 4),round(apply(estimSMSN_1[1:j,-11],2,sd),digits = 4)))
      j               <- j + 1
    }
  
  if(family == "Skew.normal")
    if(class(TrySMSN)!="try-error")
    {
      #print(TrySMSN$Sigma)
      estimSMSN_1[j,]             <- as.matrix(c(TrySMSN$beta[,1],TrySMSN$beta[,2],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$pii),byrow=TRUE)
      estimSMSN_2[(2*j-1):(2*j),] <- cbind(matrix(unlist(TrySMSN$Sigma[1]), ncol = 2, byrow = TRUE),matrix(unlist(TrySMSN$Sigma[2]), ncol = 2, byrow = TRUE))
      EP[j,]                      <- imm.smsn.mmrm(y, x, TrySMSN)$EP 
      j               <- j + 1
    }
  
}
end.timeF        <- proc.time() - start.timeI #Reset time
text             <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
########################################################################################################################################

########################################################################################################################################

########################                               Evaluacion de estimativas                                ########################

########################################################################################################################################
rbind(round(apply(estimSMSN_1,2,mean),digits = 4),t(c(betas1,betas2,shape1,shape2,pii,nu))) 

rbind(round(apply(EPSMSN_1,2,mean),digits = 4),round(apply(estimSMSN_1[,-12],2,sd),digits = 4)) 


sigmaSimul1   <- matrix(0,2,2) 
sigmaSimul2   <- matrix(0,2,2) 

for(i in seq(1,replicate*p,2))
  sigmaSimul1 <- sigmaSimul1 + estimSMSN_2[i:(i+1),1:2]
sigmaSimul1   <- sigmaSimul1/replicate; cbind(round(sigmaSimul1,digits = 4),Sigma1) 

for(i in seq(1,replicate*p,2))
  sigmaSimul2 <- sigmaSimul2 + estimSMSN_2[i:(i+1),3:4]
sigmaSimul2   <- sigmaSimul2/replicate; cbind(round(sigmaSimul2,digits = 4),Sigma2) 

round(apply(estimSMSN_1,2,sd),digits = 4)
round(apply(EPSMSN_1,2,mean),digits = 4)
########################################################################################################################################
#Manual
betas = NULL; Sigma = NULL; shape = NULL; pii = NULL; g = 2; get.init = TRUE; criteria = TRUE
group = FALSE; error = 0.0001; iter.max = 100; uni.Gama = FALSE; calc.im=FALSE
obs.prob= FALSE; kmeans.param = NULL


#########################################################################################################################################
