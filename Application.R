######################################################################
### DATA SET
#####################################################################
#Load R functions
#Ubuntu
source("densSMSN.R")
source("SMSN.mmixregre.R")
source("genmixsmsn.R")


library("EMMIXuskew") #install.packages("EMMIXuskew")
data("ais")

ais$BMI[100]   = ais$BMI[100] - mean(ais$BMI)
ais$Bfat[100]  = ais$Bfat[100] - mean(ais$Bfat)
ais$ais$Ferr[100]   = ais$ais$Ferr[100] + mean(ais$Ferr)
ais$Hg[100]  = ais$Hg[100] + mean(ais$Hg)

ais$BMI[50]   = ais$BMI[50] - mean(ais$BMI)
ais$Bfat[50]  = ais$Bfat[50] - mean(ais$Bfat)
ais$ais$Ferr[50]   = ais$ais$Ferr[50] + mean(ais$Ferr)
ais$Hg[50]  = ais$Hg[50] + mean(ais$Hg)

dataset = data.frame(ais$BMI,ais$Bfat,ais$Hg,ais$Ferr)

pairs.panels(dataset, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
pairs.panels(ais[,2:12], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

 
p   <- 2
n   <- nrow(y)
 

y  <- cbind(ais$BMI, ais$Bfat)
 
x  <- kronecker(cbind(rep(1,n),ais$Hg,ais$WCC),matrix(1,ncol(y),1))

hist(ais$BMI,main="", ylab="Density",xlab=expression(paste(y[1])), freq = FALSE) 
hist(ais$Bfat,main="", ylab="Density",xlab=expression(paste(y[2])), freq = FALSE)

respST <- smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Skew.t", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

respSN <- smsnReg.mmix(y, x, nu=3,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                     group = FALSE, family = "Skew.normal", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                     obs.prob= FALSE, kmeans.param = NULL)

respSCN<-smsnReg.mmix(y, x, nu=c(0.1,0.1),  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,
                     group = FALSE, family = "Skew.cn", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                     obs.prob= FALSE, kmeans.param = NULL)

respSSL<-smsnReg.mmix(y, x, nu=3,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,
                       group = FALSE, family = "Skew.slash", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

respN<-smsnReg.mmix(y, x, nu=3,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                     group = FALSE, family = "Normal", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                     obs.prob= FALSE, kmeans.param = NULL)

respT<-smsnReg.mmix(y, x, nu=3,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,
                    group = FALSE, family = "t", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                    obs.prob= FALSE, kmeans.param = NULL)


round(respSN$beta, digits=4)
respSN$Sigma
respSN$shape
respSN$nu
respSN$pii


plot(x1[,2],y[,1], col = (respT$group==1),ylab="Math",xlab="Writing", pch=16)
par(new = TRUE)
plot(x1[,2],y[,1], col = (respT$group==2),ylab="Math",xlab="Writing",pch = 1)
abline(respT$beta[[1]][1], respT$beta[[1]][2])
abline(respT$beta[[2]][1], respT$beta[[2]][2])
abline(respN$beta[[1]][1], respN$beta[[1]][2],lty=2)
abline(respN$beta[[2]][1], respN$beta[[2]][2],lty=2)

#######################################################################
#install.packages("psych")
library(psych)

dat.csv              <- read.csv("hsb2-2.csv")



head(dat.csv)
dat.csv$math[100]    = dat.csv$math[100] - 40
dat.csv$science[100] = dat.csv$science[100] - 40
dat.csv$write[100]   = dat.csv$write[100] + mean(dat.csv$write)
dat.csv$read[100]    = dat.csv$read[100] - mean(dat.csv$read)

dat.csv$math[150]    = dat.csv$math[150] - 50
dat.csv$science[150] = dat.csv$science[150] - 50
dat.csv$write[150]   = dat.csv$write[15] + mean(dat.csv$write)
dat.csv$read[150]    = dat.csv$read[15] - mean(dat.csv$read)

########################################################3

dat.csv$math[100]    = dat.csv$math[100] - mean(dat.csv$math)
dat.csv$science[100] = dat.csv$science[100] - mean(dat.csv$science)
dat.csv$write[100]   = dat.csv$write[100] + mean(dat.csv$write)
dat.csv$read[100]    = dat.csv$read[100] + mean(dat.csv$read)
dat.csv$socst[100]   = dat.csv$socst[100] + mean(dat.csv$socst)


dat.csv$math[150]    = dat.csv$math[150] - mean(dat.csv$math)
dat.csv$science[150] = dat.csv$science[150] - mean(dat.csv$science)
dat.csv$write[150]   = dat.csv$write[150] + mean(dat.csv$write)
dat.csv$read[150]    = dat.csv$read[150] + mean(dat.csv$read)
dat.csv$socst[150]   = dat.csv$socst[150] + mean(dat.csv$socst)



pairs.panels(dat.csv[,7:11], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pairs.panels(dat.csvSO[,7:11], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

y  <- cbind( dat.csv$math,dat.csv$science)
#x1 <- cbind(1,dat.csv$write)
#x2 <- cbind(1,dat.csv$read)
 
x  <- list(kronecker(x1,matrix(1,ncol(y),1)),kronecker(x2,matrix(1,ncol(y),1)))

###############################################################################################

fitST   <- smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Skew.t", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

fitST_2 <- smsnReg.mmix(y_2, x_2, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                      group = FALSE, family = "Skew.t", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                      obs.prob= FALSE, kmeans.param = NULL)

fitSN   <- smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Skew.normal", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

fitSN_2   <- smsnReg.mmix(y_2, x_2, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                        group = FALSE, family = "Skew.normal", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                        obs.prob= FALSE, kmeans.param = NULL)


fitSCN <- smsnReg.mmix(y, x, nu=c(0.1,0.1),  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                     group = FALSE, family = "Skew.cn", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                     obs.prob= FALSE, kmeans.param = NULL)

fitSCN_2 <- smsnReg.mmix(y_2, x_2, nu=c(0.1,0.1),  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Skew.cn", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

fitSSL <- smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Skew.slash", error = 10^-5, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

fitSSL_2 <- smsnReg.mmix(y_2, x_2, nu=3,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Skew.slash", error = 10^-4, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

fitN <- smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,
                       group = FALSE, family = "Normal", error = 10^-5, iter.max = 200, uni.Gama = FALSE, calc.im=FALSE,
                       obs.prob= FALSE, kmeans.param = NULL)

fitT <- smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,
                     group = FALSE, family = "t", error = 10^-5, iter.max = 200, uni.Gama = FALSE, calc.im=FALSE,
                     obs.prob= FALSE, kmeans.param = NULL)

 
