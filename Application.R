######################################################################
### DATA SET
#####################################################################
#Load R functions
#Ubuntu
source("densSMSN.R")
source("SMSN.mmixregre.R")
source("genmixsmsn.R")


#Dados
dat.csv              <- read.csv("hsb2-2.csv")

y  <- cbind( dat.csv$science,dat.csv$write)
x  <- kronecker(cbind(rep(1,n),dat.csv$read),matrix(1,ncol(y),1))

#Figure 6
plot(density(y[,1]),main="", ylab="Density",xlab=expression(paste(y[1]))) #a)
plot(density(y[,2]),main="", ylab="Density",xlab=expression(paste(y[2]))) #b)
plot(density(dat.csv$read),main="", ylab="Density",xlab=expression(paste(x[1]))) #c)


###############################################################################################

initial.values.SN  <- initial.Values.fm.smsn.mr(y,x,g=2,get.init="k-means",family="Skew.normal",lower=1.01,upper=20,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
fitSN_g2           <- smsnReg.mmix(y, x, nu=4,  betas = initial.values.SN$beta, Sigma = initial.values.SN$Sigma, shape = initial.values.SN$shape, pii = initial.values.SN$pii, g = 2, get.init = FALSE, criteria = TRUE,group = FALSE, family = "Skew.normal", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,obs.prob= FALSE, kmeans.param = NULL)

initial.values.ST  <- initial.Values.fm.smsn.mr(y,x,g=2,get.init="k-means",family="Skew.t",lower=1.01,upper=20,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
fitST_g2           <- smsnReg.mmix(y, x, nu=initial.values.ST$nu,  betas = initial.values.ST$beta, Sigma = initial.values.ST$Sigma, shape = initial.values.ST$shape, pii = initial.values.ST$pii, g = 2, get.init = FALSE, criteria = TRUE,group = FALSE, family = "Skew.t", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,obs.prob= FALSE, kmeans.param = NULL)

initial.values.SCN <- initial.Values.fm.smsn.mr(y,x,g=2,get.init="k-means",family="Skew.cn",lower=0.1,upper=0.9,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
fitSCN_g2          <- smsnReg.mmix(y, x, nu=initial.values.SCN$nu,  betas = initial.values.SCN$beta, Sigma = initial.values.SCN$Sigma, shape = initial.values.SCN$shape, pii = initial.values.SCN$pii, g = 2, get.init = FALSE, criteria = TRUE,group = FALSE, family = "Skew.cn", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,obs.prob= FALSE, kmeans.param = NULL)

EP_SN  = imm.smsn.mmrm(y, x, fitSN_g2)$EP
EP_ST  = imm.smsn.mmrm(y, x, fitST_g2)$EP
EP_SCN = imm.smsn.mmrm(y, x, fitSCN_g2)$EP

 
