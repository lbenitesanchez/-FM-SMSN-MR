#Load packages
library("magic")#install.packages("magic")
library("mixsmsn")
library(mnormt) #install.packages("mnormt")

 
#Configuration
replicate   <- 1
p           <- 2
family      <- "Skew.t"


if(family == "Skew.t" || family == "Skew.slash")
{
 ## Configuracion de objetos para familia Skew.t o Skew.slash
 estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3 + 2 + 2*2 + 2 + 1)))
 EPSMSN_1              <- array(0,dim=c(replicate,ncol=(3 + 2 + 2*2 + 2)))
 colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2","nu")
 colnames(EPSMSN_1)    <- c("EPBeta1.1","EPBeta1.2","EPBeta1.3","EPBeta2.1","EPBeta2.2","EPShape1.1","EPShape1.2","EPShape2.1","EpShape2.2","EPpii1","EPnu")
 estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
 EPSMSN_2              <- array(0,dim=c(replicate*p,ncol=4))
}

if(family == "Skew.cn")
{
 ## Configuracion de objetos para familia Skew.t o Skew.slash
 estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3 + 2 + 2*2 + 2 + 1 + 1)))
 colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2","nu","gama")
 estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
}

if(family == "Skew.normal")
{
 ## Configuracion de objetos para familia Skew.normal (no tiene grados de libertad)
 estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3*2 + 2*2 + 1)))
 colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2")
 estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
}

##################################################################################################################
## Gerando dados mixturas de regressoes smsn multi
##################################################################################################################

pii     <- c(0.3,0.7)
n       <- 1000 #Sample size
x1      <- cbind(rep(1,n*p),runif(n*p,0,1),runif(n*p, 2,4))
x2      <- cbind(rep(1,n*p),runif(n*p,0,1))
x       <- list(x1,x2)

#Freedom degrees
if(family == "Skew.t" || family == "Skew.slash")
{ nu   <- 3   #No rodar para Skew.normal e normal
}else{
  nu   <- c(0.1,0.1) #Para Skew.cn
}


if(family == "Skew.t")      {k1  <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2);b <- -sqrt(2/pi)*k1}
if(family == "Skew.normal") {k1  <- 1                                     ;b <- -sqrt(2/pi)*k1}
if(family == "Skew.cn")     {k1  <- nu[1]/nu[2]^(1/2)+1-nu[1]             ;b <- -sqrt(2/pi)*k1}
if(family == "Skew.slash")  {k1  <- 2*nu/(2*nu-1)                         ;b <- -sqrt(2/pi)*k1}

#First component
 betas1 <- c(-1,-4,-3) 
 Sigma1 <- matrix(c(3,1,1,3), 2,2)
 shape1 <- matrix(c(3,5),2,1)
 delta1 <- shape1/sqrt(1 + sum(shape1^2))
 Delta1 <- matrix.sqrt(Sigma1)%*%delta1
 mu1    <- t(matrix(x[[1]]%*%betas1 + rep(b*Delta1,n),p,n))

#Second component
 betas2 <- c(3,7)   
 Sigma2 <- matrix(c(2,1,1,2), 2,2)
 shape2 <- matrix(c(1,4),2,1)
 delta2 <- shape2/(sqrt(1 + sum(shape2^2)))
 Delta2 <- matrix.sqrt(Sigma2)%*%delta2
 mu2    <- t(matrix(x[[2]]%*%betas2+rep(b*Delta2,n),p,n))

########################################################################################################################################

 if(family == "Skew.t" || family == "Skew.slash" || family == "Skew.cn")
 {     y <- genmixsmsn(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu);hist(y[,1])
 }else{y <- genmixsmsn(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu=NULL);hist(y)}

 TrySMSN      <- try(smsnReg.mmix(y, x, nu,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,group = FALSE, family = family, error = 10^-6, iter.max = 100, uni.Gama = FALSE, calc.im=FALSE, obs.prob= FALSE, kmeans.param = NULL))
 estimSMSN_1  <- as.matrix(c(TrySMSN$beta[[1]],TrySMSN$beta[[2]],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$pii,TrySMSN$nu),byrow=TRUE)
 
 G   <- c(sum(TrySMSN$group==1),sum(TrySMSN$group==2))

 #Load package
 library(ggplot2)
 library(lattice)
 
dat        <- cbind(y,x[[1]][1:1000,2:3],as.factor(c(rep(1,G[1]),rep(2,G[2]))))
dat        <- data.frame(dat)
names(dat) <- c("y1","y2","x1","x2","Group")
dat$Group  <- as.factor(dat$Group)
 
postscript("clusteringScenario1.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
 ggplot(dat, aes(y1, y2, Group, color = Group)) + geom_point() + stat_ellipse() +
   theme_bw() +
   theme(axis.text = element_text(size = 14),
         #legend.key = element_rect(fill = "navy"),
         legend.background = element_rect(fill = "white"),
         legend.position = c(0.9, 0.8),
         #panel.grid.major = element_line(colour = "grey40"),
         panel.grid.minor = element_blank()
        )
dev.off() #Fechando o dispositivo potscript

postscript("clusteringScenario1_group2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
 xyplot(cbind(y1,y2) ~ x[[1]][1:1000,2:3], dat, grid = TRUE, group = Group, pch=c(19,16),xlab="", ylab="", main="Real Data" ,col=c("red","cyan"))  
dev.off
######
dat2        <- cbind(y,x[[2]][1:1000,2],as.factor(c(rep(1,G[1]),rep(2,G[2]))))
dat2        <- data.frame(dat2)
names(dat2) <- c("y1","y2","x1","Group")
dat2$Group   <- as.factor(dat2$Group)

postscript("clusteringScenario1_group2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
ggplot(dat2, aes(y1, y2, Group, color = Group)) + geom_point() + stat_ellipse() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        #legend.key = element_rect(fill = "navy"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.9, 0.8),
        #panel.grid.major = element_line(colour = "grey40"),
        panel.grid.minor = element_blank()
  )
dev.off() #Fechando o dispositivo potscript

postscript("clusteringScenario1_group2.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
xyplot(cbind(y1,y2) ~ x[[2]][1:1000,2], dat, grid = TRUE, group = Group, pch=c(19,16), xlab="", ylab="", main="Real Data" ,col=c("red","cyan"))  
dev.off
########################################################################################################################################
