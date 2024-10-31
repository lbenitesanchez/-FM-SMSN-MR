 
#Load packages
library("magic")#install.packages("magic")
library("mixsmsn")
library(mnormt) #install.packages("mnormt")

#Load R functions

#Configuration
replicate   <- 500
p           <- 2
family      <- "Skew.t"


if(family == "Skew.t" || family == "Skew.slash")
{
  ## Configuracion de objetos para familia Skew.t o Skew.slash
  estimSMSN_1           <- array(0,dim=c(replicate,ncol=(3 + 2 + 2*2 + 2 + 1)))
  vies                  <- array(0,dim=c(replicate,ncol=(3 + 2 + 2*2 + 2)))
  rmse                  <- array(0,dim=c(replicate,ncol=(3 + 2 + 2*2 + 2)))
  colnames(estimSMSN_1) <- c("Beta1.1","Beta1.2","Beta1.3","Beta2.1","Beta2.2","Shape1.1","Shape1.2","Shape2.1","Shape2.2","pii1","pii2","nu")
  colnames(vies)        <- c("ViesBeta1.1","ViesBeta1.2","ViesBeta1.3","ViesBeta2.1","ViesBeta2.2","ViesShape1.1","ViesShape1.2","ViesPShape2.1","ViesShape2.2","Viespii1","Viesnu")
  colnames(vies)        <- c("rmseBeta1.1","rmseBeta1.2","rmseBeta1.3","rmseBeta2.1","rmseBeta2.2","rmseShape1.1","rmseShape1.2","rmsePShape2.1","rmseShape2.2","rmsepii1","rmsenu")
  estimSMSN_2           <- array(0,dim=c(replicate*p,ncol=4))
  vies_2                <- array(0,dim=c(replicate*p,ncol=4))
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

vt     <- c(betas1,betas2,shape1,shape2,pii,nu)
########################################################################################################################################
j       <- 1
start.timeI      <- proc.time(); #Begin time
while(j <= replicate)
{
  cat("#Replicate",j,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')
  
  if(family == "Skew.t" || family == "Skew.slash" || family == "Skew.cn")
  {     y <- genmixsmsn(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu);hist(y[,1])
  }else{y <- genmixsmsn(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu=NULL);hist(y)}
  
  TrySMSN <- try(smsnReg.mmix(y, x, nu,  betas = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = 2, get.init = TRUE, criteria = TRUE,group = FALSE, family = family, error = 10^-6, iter.max = 500, uni.Gama = FALSE, calc.im=FALSE, obs.prob= FALSE, kmeans.param = NULL))
  
  if(family == "Skew.t" || family == "Skew.slash" || family == "Skew.cn")
    if(class(TrySMSN)!="try-error" && TrySMSN$beta[[1]][1]<0 && TrySMSN$beta[[1]][1] > -3 &&  unlist(TrySMSN$shape[2])[1]>0 &&  unlist(TrySMSN$shape[2])[2]>0  &&  unlist(TrySMSN$shape[1])[1]<4)
    {
      estimSMSN_1[j,]             <- as.matrix(c(TrySMSN$beta[[1]],TrySMSN$beta[[2]],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$pii,TrySMSN$nu),byrow=TRUE)
      estimSMSN_2[(2*j-1):(2*j),] <- cbind(matrix(unlist(TrySMSN$Sigma[1]), ncol = 2, byrow = TRUE),matrix(unlist(TrySMSN$Sigma[2]), ncol = 2, byrow = TRUE))
      j               <- j + 1
    }
  
  if(family == "Skew.normal")
    if(class(TrySMSN)!="try-error" && TrySMSN$beta[[1]][1]<0)
    {
      estimSMSN_1[j,]             <- as.matrix(c(TrySMSN$beta[[1]],TrySMSN$beta[[2]],unlist(TrySMSN$shape[1]),unlist(TrySMSN$shape[2]),TrySMSN$pii),byrow=TRUE)
      estimSMSN_2[(2*j-1):(2*j),] <- cbind(matrix(unlist(TrySMSN$Sigma[1]), ncol = 2, byrow = TRUE),matrix(unlist(TrySMSN$Sigma[2]), ncol = 2, byrow = TRUE))
      j               <- j + 1
    }
}

vies             <- (1/replicate)*(apply(estimSMSN_1-matrix(vt,nr=replicate,nc=ncol(estimSMSN_1),byrow=T),2,sum))
rmse             <- sqrt(apply((estimSMSN_1-matrix(vt,nr=replicate,nc=ncol(estimSMSN_1),byrow=T))^2,2,mean)) #vt: valores teoricos

Sigma            <- cbind(Sigma1,Sigma2)
SigmaSuma1       <- estimSMSN_2-do.call("rbind", rep(list(Sigma), replicate))
SigmaSuma2       <- (estimSMSN_2-do.call("rbind", rep(list(Sigma), replicate)))^2
SigmaSuma2_1     <- 0
SigmaSuma2_2     <- 0
for(j in 1:replicate)
{
  SigmaSuma2_1 <- SigmaSuma2_1 + SigmaSuma1[(2*j-1):(2*j),]
  SigmaSuma2_2 <- SigmaSuma2_2 + SigmaSuma2[(2*j-1):(2*j),] 
}

vies2            <- (1/replicate)*SigmaSuma2_1
rmse2            <- sqrt(SigmaSuma2_2/replicate)

end.timeF        <- proc.time() - start.timeI #Reset time
text             <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")


#####################################################################################################################3
#                                                    Graficos
#####################################################################################################################3
gray1 <- rgb(51,51,51,max=255)
gray2 <- rgb(112,112,112,max=255)
gray3 <- rgb(150,150,150,max=255)

n <- c(100,250,500,100,2500)
#Vies  

VIES_SN = rbind(c(100, -0.36717311 , 0.03383585 , 0.11338743 ,   0.02354057, -0.02123391, -1.79043019, -1.10980917, -0.84252660, -1.07405314, 0.00162000),
                c(250, -0.077532286, 0.011757031, 0.024946306,  0.004819719, -0.005853751, -1.905473108, -2.138022035, -0.832592492, -1.779127638,  0.000424000),
                c(500, -0.064620873, 0.013706654, 0.020984310, -0.001780419, -0.004011008, -2.045661487, -3.214992656, -0.349537201, -0.970509466, -0.000280000 ),
                c(1000,-0.005929143, 0.007187810, 0.001269273,  0.004680149, -0.007473595, -1.547945625, -2.192136082, -0.536280883, -2.108766270, -0.000136000),
                c(2500,     0.00328,      0.0039,      0.0017,      0.00400,  -0.0034,   0.0435,   0.0504,   0.1048,   0.2024,   0.0002)) 

VIES_SN_Sigma = rbind(c(100,  0.3274983,  2.3274983,  -3.286205,  0.927498,  -1.086205,  -1.286205),
                      c(250, -0.03586547, 1.96413453, -2.8765104, 0.7641345, -0.6765104, -0.8765104),
                      c(500,  -0.04623817,1.95376183, -2.7359822, 0.6537618, -0.5359822, -0.7359822),
                      c(1000, -0.0795102, 1.9204898, -2.7116666, 0.5204898, -0.5116666, -0.7116666),
                      c(2500, -0.1330387, 1.8669613, -2.5073295, 0.3669613, -0.3073295, -0.5073295)) 

VIES_ST = rbind(c(100, -0.356879000,  0.015179570,  0.073192913,  0.054941774,  0.011839583, -2.745318829,  2.713496702,  1.020372974,  2.967551490,  0.002399724, 0.610229519),
                c(250, -0.082863060, -0.002115239,  0.015069568,  0.045180080,  0.002642348, -0.669530599, -0.497413129,  0.204770402,  0.526624247,  0.002040765, 0.189415398),
                c(500, -0.004789173, -0.001625126,  0.004523428,  0.02120762, -0.000379370, -0.289152310, -0.285334354,  0.148087821,  0.267811873,  0.000719784, 0.113274638),
                c(1000, 0.050046505, -0.002823283, -0.005817362,  0.028434102, -0.005806229, -0.121626522, -0.138414382,  0.104630528,  0.199593793,  0.000207596, 0.126600251),
                c(2500,      0.0328,       0.0039,       0.0017,   0.00400,  -0.0034,   0.0435,   0.0504,   0.1048,   0.2024,   0.0002,   0.0981)) 

VIES_ST_Sigma = rbind(c(100,   -0.3486788, -0.8692656,  -0.5582268,   0.2451514,  0.1611563,  0.1154942),
                      c(250,   -0.1753618, -0.2931099,  -0.1916233,   0.1180455,  0.07919030, 0.09454698),
                      c(500,  0.002959826, -0.100894251, 0.006253244, 0.1174482,  0.0954346,  0.1204875),
                      c(1000,  0.10107580, 0.02514546,   0.14509763,  0.1128477,  0.1054036,  0.1414460),
                      c(2500,  0.1538978,  0.1084130,    0.2008392,   0.09182983, 0.08921372, 0.12865850)) 

VIES_SCN= rbind(c(100, -0.3832571778,  0.0399156746,  0.0741980724, -0.0296468205,  0.0127115089, -2.7926604468, -0.8945936738,  0.1608639453, -0.1672250970,  0.0008550275, 0.0138194542, -0.0004944827),
                c(250, -0.0925675278,  0.0334931111,  0.0108211419, -0.0048198639,  0.0019156780, -0.9928197243, -0.9757741807,  0.0558965321, -0.0699991748, -0.0004290252, 0.0038413123, 0.0012962385),
                c(500, -0.0193780378, -0.0139445312,  0.0047286366, -0.0035410082, -0.0060714606, -0.4320466393, -0.4646150031, -0.0285164098, -0.1533873652, -0.0009212453, 0.0011957100, 0.0026980284),
                c(1000, 8.040137e-03, -6.433969e-03, -1.852360e-03, -5.011103e-04, -3.880585e-04, -2.686117e-01, -4.460724e-01, -2.586989e-02, -1.980702e-01, -4.356193e-04, 8.979032e-04, -7.544087e-06),
                c(2500, 0.0051519161,  0.0013867308, -0.0011722187, -0.0016845897,  0.0011779741, -0.2149680060, -0.4051837185, -0.0336267947, -0.1909875114, -0.0002335222, 0.0003107089, 0.0003231051)) 


VIES_SCN_Sigma = rbind(c(100,   -0.09465773, -1.1623107,  -0.8731328,  0.07113181,  -0.02340447,  -0.19161083),
                       c(250,   -0.2625262,  -0.4356594,  -0.3478464,  0.034978555, -0.004206477, -0.069590962),
                       c(500,   -0.1291759,  -0.191370,   -0.125392,   0.001029754, -0.03273435,  -0.07728641),
                       c(1000,  -0.03767504, -0.06529411, -0.09562497, 0.007674927, -0.009390895, -0.043245466),
                       c(2500,  0.003132237, -0.03410554, -0.06193580, 0.001821698, -0.009587574, -0.030390204)) 


colnames(VIES_SN)        <- c("Sample Size","ViesBeta1.1","ViesBeta1.2","ViesBeta1.3","ViesBeta2.1","ViesBeta2.2","ViesShape1.1","ViesShape1.2","ViesPShape2.1","ViesShape2.2","Viespii1")
colnames(VIES_SN_Sigma)  <- c("Sample Size","ViesSigma1.11","ViesSigma1.12","ViesSigma1.22","ViesSigma2.11","ViesSigma2.12","ViesSigma2.22")
colnames(VIES_ST)        <- c("Sample Size","ViesBeta1.1","ViesBeta1.2","ViesBeta1.3","ViesBeta2.1","ViesBeta2.2","ViesShape1.1","ViesShape1.2","ViesPShape2.1","ViesShape2.2","Viespii1","Viesnu")
colnames(VIES_ST_Sigma)  <- c("Sample Size","ViesSigma1.11","ViesSigma1.12","ViesSigma1.22","ViesSigma2.11","ViesSigma2.12","ViesSigma2.22")
colnames(VIES_SCN)       <- c("Sample Size","ViesBeta1.1","ViesBeta1.2","ViesBeta1.3","ViesBeta2.1","ViesBeta2.2","ViesShape1.1","ViesShape1.2","ViesPShape2.1","ViesShape2.2","Viespii1","Viesnu","Viesgamma")
colnames(VIES_SCN_Sigma) <- c("Sample Size","ViesSigma1.11","ViesSigma1.12","ViesSigma1.22","ViesSigma2.11","ViesSigma2.12","ViesSigma2.22")

RMSE_SN = rbind(c(100,  0.90115108, 0.60578887, 0.30949254, 0.15906567, 0.26245969, 6.01728976, 7.26969054, 2.07732125, 3.01558461, 0.04433734),
                c(250,  0.50428481, 0.30674486, 0.12636034, 0.10183405, 0.15697636, 3.27633313, 4.35758932, 1.36341841, 2.51018098, 0.03096889 ),
                c(500,  0.40414372, 0.20384503, 0.10391237, 0.06338919, 0.10386322, 2.29846846, 4.00204627, 0.70946657, 1.42524536, 0.02123356),
                c(1000, 0.20296592, 0.10909807, 0.00539035, 0.02446866, 0.07112004, 2.10770397, 3.13365290, 0.73313922, 2.57260665, 0.01438471 ),
                c(2500, 0.1436,      0.1059,     0.0333,     0.0451,     0.0466,     0.4963,     0.7151,      0.2303,      0.4409,     0.0091)) 

RMSE_SN_Sigma = rbind(c(100,  1.490308,  2.744269,   3.418037, 1.968760, 1.593160, 1.471957),
                      c(250,  0.8321237, 2.1328310,  2.980325, 1.273067, 1.173157, 1.032072),
                      c(500,  0.6076402, 2.0455498,  2.836958, 1.129934, 1.050905, 1.591128),
                      c(1000, 0.4456196, 1.9699076, 2.8005455, 1.019586, 0.9981927, 0.8492490),
                      c(2500, 0.3240123, 1.8901929, 2.5846868, 0.9159184, 0.8070241, 0.6327728))  

RMSE_ST = rbind(c(100,  0.81271106,  0.73214888, 0.27242739, 0.27209690, 0.29355749, 5.37639024, 11.32596251, 2.89385858,  8.25088421, 0.04625682, 4.68289421),
                c(250,  0.55354992,  0.36102005, 0.16418578, 0.17321480, 0.19225656, 1.13811672, 1.80286836,  0.71194399,  1.57464565, 0.02822353, 0.59477539 ),
                c(500,  0.39075604,  0.24932275, 0.11910583, 0.12868264, 0.12901100, 0.78410985, 1.25881356,  0.48628767,  0.97553109, 0.02005792, 0.37936359),
                c(1000, 0.30368799,  0.17336888, 0.08555580, 0.09473368, 0.08844196, 0.63481754, 0.94457633,  0.34230507,  0.67580062, 0.01448885, 0.30097110 ),
                c(2500, 0.1836,      0.1059,     0.0533,     0.0651,     0.0566,     0.4963,     0.7151,      0.2303,      0.4709,     0.0091,     0.1906 )) 

RMSE_ST_Sigma = rbind(c(100,  1.209997,   1.219282,  1.282156, 0.7805734, 0.6073919, 0.6821505),
                      c(250,  0.8346481, 0.5957959, 0.8257699, 0.4596543, 0.3848933, 0.4566116),
                      c(500,  0.6054077, 0.4155354, 0.5765562, 0.3441158, 0.3067591, 0.3567213),
                      c(1000, 0.4703202, 0.3305961, 0.4968067, 0.2522892, 0.2202212, 0.2719938),
                      c(2500, 0.3356564, 0.2409532, 0.3775837, 0.1699129, 0.1586959, 0.1991748)) 

RMSE_SCN= rbind(c(100,  0.84554448, 0.70318485,   0.25887929, 0.21015002, 0.28565544, 4.13280397, 5.59213149, 0.93062620, 1.67881160, 0.04426347, 0.07989008, 0.04665162),
                c(250,  0.52743105, 0.33966758,   0.15499434, 0.11821421, 0.16963197, 1.58234368, 2.14589348, 0.53143490, 1.02397699, 0.02905290, 0.03353628, 0.02741510 ),
                c(500,  0.38830565, 0.21150455,   0.12153861, 0.08095662, 0.11971932, 0.86554723, 1.41492260, 0.38976263, 0.76768867, 0.01890508, 0.02230516, 0.01864529),
                c(1000, 0.25301040, 0.16041719,   0.07895480, 0.05904225, 0.08581814, 0.65365350, 1.04535115, 0.26582662, 0.53183185, 0.01429892, 0.01536692, 0.01325814 ),
                c(2500, 0.161249018, 0.098418572, 0.050016447, 0.036809938, 0.052361325, 0.480566469, 0.720163060, 0.173937237, 0.375194353, 0.009595265, 0.008541815, 0.007359936)) 

RMSE_SCN_Sigma = rbind(c(100,  1.478835, 1.468523,  1.393027,  0.6232261,  0.4574647, 0.5537762),
                       c(250,  0.7433152, 0.7487173, 0.8067583, 0.3568948, 0.2971856, 0.3084168),
                       c(500,  0.5480918, 0.3939318, 0.5503424, 0.2358940, 0.2120977, 0.2346825),
                       c(1000, 0.3953824, 0.2722044, 0.3857456, 0.1707818, 0.1587042, 0.1581789),
                       c(2500, 0.2506759, 0.1670999, 0.2307841, 0.1114429, 0.1024816, 0.1069983)) 

colnames(RMSE_SN)        <- c("Sample Size","rmseBeta1.1","rmseBeta1.2","rmseBeta1.3","rmseBeta2.1","rmseBeta2.2","rmseShape1.1","rmseShape1.2","rmsePShape2.1","rmseShape2.2","rmsepii1")
colnames(RMSE_SN_Sigma)  <- c("Sample Size","RmsesSigma1.11","RmseSigma1.12","RmseSigma1.22","RmseSigma2.11","RmseSigma2.12","RmseSigma2.22")
colnames(RMSE_ST)        <- c("Sample Size","rmseBeta1.1","rmseBeta1.2","rmseBeta1.3","rmseBeta2.1","rmseBeta2.2","rmseShape1.1","rmseShape1.2","rmsePShape2.1","rmseShape2.2","rmsepii1","rmsenu")
colnames(RMSE_ST_Sigma)  <- c("Sample Size","RmsesSigma1.11","RmseSigma1.12","RmseSigma1.22","RmseSigma2.11","RmseSigma2.12","RmseSigma2.22")
colnames(RMSE_SCN)       <- c("Sample Size","rmseBeta1.1","rmseBeta1.2","rmseBeta1.3","rmseBeta2.1","rmseBeta2.2","rmseShape1.1","rmseShape1.2","rmsePShape2.1","rmseShape2.2","rmsepii1","rmsenu","rmsegamma")
colnames(RMSE_SCN_Sigma) <- c("Sample Size","RmsesSigma1.11","RmseSigma1.12","RmseSigma1.22","RmseSigma2.11","RmseSigma2.12","RmseSigma2.22")

 