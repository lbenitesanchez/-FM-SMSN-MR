Nboots<-500  ## Numero de amostras Bootstrap
#################################
 

#################################
# Normal
family      <- "Skew.normal"

resB  <- array(0,dim=c(Nboots,ncol=9))
resB2 <- array(0,dim=c(Nboots,ncol=6))
j     <- 1
while (j <= Nboots){
  print(j)
  pp      <- sample(seq(1,n), n, replace = TRUE)
  y_1     <- y[pp,]
  x_1     <- x[  c(rbind(pp, pp + 1)),]
  initial.values.SN  <- initial.Values.fm.smsn.mr(y_1,x_1,g=2,get.init="k-means",family="Skew.normal",lower=1.01,upper=20,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
  TrySMSN           <- try(smsnReg.mmix(y_1, x_1, nu=4,  betas = initial.values.SN$beta, Sigma = initial.values.SN$Sigma, shape = initial.values.SN$shape, pii = initial.values.SN$pii, g = 2, get.init = FALSE, criteria = TRUE,group = FALSE, family = "Skew.normal", error = 10^-8, iter.max = 300, uni.Gama = FALSE, calc.im=FALSE,obs.prob= FALSE, kmeans.param = NULL))
  
  
  if(class(TrySMSN)!="try-error")
    if(TrySMSN$pii[1]<0.5 && abs(TrySMSN$shape[[1]][2]) < 10 && abs(TrySMSN$shape[[2]][2]) < 10 )
    {
      resB[j,]  <- c(TrySMSN$beta[,1],TrySMSN$beta[,2],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$pii[1])                
      resB2[j,] <- c((unlist(TrySMSN$Sigma[1]))[1],(unlist(TrySMSN$Sigma[1]))[2],(unlist(TrySMSN$Sigma[1]))[4],(unlist(TrySMSN$Sigma[2]))[1],(unlist(TrySMSN$Sigma[2]))[2],(unlist(TrySMSN$Sigma[2]))[4] )
      j         <- j + 1
      #print(TrySMSN$beta)
    }
}                                                   

EP.B <- round(apply(resB,2,sd),digits=4)
EP.B

EP.B2 <- round(apply(resB2,2,sd),digits=4)
EP.B2

#################################
family      <- "Skew.t"
resB  <- array(0,dim=c(Nboots,ncol=12))
resB2 <- array(0,dim=c(Nboots,ncol=6))
j     <- 1
while (j <= Nboots){
  print(j)
  pp   <- sample(seq(1,n), n, replace = TRUE)
  y_1   <- y[pp]
  x_1   <- x[pp,]
  TrySMSN <- try(smsnReg.mmix(y, x, nu=4,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,group = FALSE, family = family, error = 10^-6, iter.max = 100, uni.Gama = FALSE, calc.im=FALSE, obs.prob= FALSE, kmeans.param = NULL))
  if(class(TrySMSN)!="try-error" && TrySMSN$pii[1]>0.9)
  {
    resB[j,]  <- c(TrySMSN$beta[,1],TrySMSN$beta[,2],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$nu,TrySMSN$pii[1])                
    resB2[j,] <- c((unlist(TrySMSN$Sigma[1]))[1],(unlist(TrySMSN$Sigma[1]))[2],(unlist(TrySMSN$Sigma[1]))[4],(unlist(TrySMSN$Sigma[2]))[1],(unlist(TrySMSN$Sigma[2]))[2],(unlist(TrySMSN$Sigma[2]))[4] )
    j         <- j + 1
  }
}                                                   

EP.B <- round(apply(resB,2,sd),digits=4)
EP.B


##################################################################
Nboots <- 500  ## Numero de amostras Bootstrap
family <- "Skew.cn"
resB   <- array(0,dim=c(Nboots,ncol=13))
resB2  <- array(0,dim=c(Nboots,ncol=6))
j      <- 1
while (j <= Nboots){
  print(j)
  pp   <- sample(seq(1,n), n, replace = TRUE)
  y_1   <- y[pp,]
  x_1   <- x[c(rbind(pp, pp + 1)),]
  TrySMSN <- try(smsnReg.mmix(y_1, x_1, nu=c(0.1,0.1),  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = FALSE,group = FALSE, family = family, error = 10^-5, iter.max = 200, uni.Gama = FALSE, calc.im=FALSE, obs.prob= FALSE, kmeans.param = NULL))
  if(class(TrySMSN)!="try-error")
  if(TrySMSN$pii[1]>0.8)
  {
    resB[j,]  <- c(TrySMSN$beta[,1],TrySMSN$beta[,2],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$nu,TrySMSN$pii[1])                
    resB2[j,] <- c((unlist(TrySMSN$Sigma[1]))[1],(unlist(TrySMSN$Sigma[1]))[2],(unlist(TrySMSN$Sigma[1]))[4],(unlist(TrySMSN$Sigma[2]))[1],(unlist(TrySMSN$Sigma[2]))[2],(unlist(TrySMSN$Sigma[2]))[4] )
    j         <- j + 1
  }
}                                                   

EP.B <- round(apply(resB,2,sd),digits=4)
EP.B

EP.B2 <- round(apply(resB2,2,sd),digits=4)
EP.B2
