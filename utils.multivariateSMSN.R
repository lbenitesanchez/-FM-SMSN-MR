initial.Values.fm.smsn.mr <- function(y,x,g=2,get.init="k-means",family="Skew.t",lower=1,upper=20,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
{
  p  <- ncol(y)
  n  <- nrow(y)
  q  <- ncol(x)
  
  ######################################################################################
  #Segue os tres tipos de algoritmos de particao que podem ser utilizados para obter
  #os chutes iniciais para beta_j, sigma_j e pii_j
  
  
  if(get.init == "k-means")
  {
    if(length(g) == 0) stop("g is not specified correctly.\n")
    
    k.iter.max     <- 200
    n.start        <- 1
    algorithm      <- "Hartigan-Wong"
    
    if(g > 1)
    {
      init         <- kmeans(y,g,k.iter.max,n.start,algorithm)
      pii          <- init$size/nrow(y)
      shape        <- Sigma <- list()
      betas        <- matrix(0,q,g)
      cc           <- init$cluster
      for (j in 1:g)
      {
        betas[,j]            <- solve(t(x[cc==j,])%*%x[cc==j,])%*%t(x[cc==j,])%*%c(t(y[cc==j])) #init$centers[j,]
        shape[[j]]           <- 3*sign(apply((y[cc == j, ] - t(matrix(t(x[cc==j,]%*%betas[,j]),p, nrow(y[cc == j, ]))))^3, 2, sum))
        Sigma[[j]]           <- var(y[cc == j,])
        dimnames(Sigma[[j]]) <- NULL
        names(shape[[j]])    <- NULL
      }
    }else{
      pii                    <- 1
      betas                  <- solve(t(x)%*%x)%*%t(x)%*%c(t(y))
      Sigma[[1]]             <- var(y)
      shape[[1]]             <- 3*sign(apply((y - t(matrix(t(x%*%betas),p, nrow(y))))^3, 2, sum))
      dimnames(Sigma[[1]])   <- NULL
      names(mu[[1]])         <- NULL
      names(shape[[1]])      <- NULL
    }
    initial_values <- list(beta = betas, Sigma = Sigma, shape = shape, pii = pii)#;print(initial_values)
  }
  
 
  ########################################################################################################################
  #-----------------------------------------------------------------------------------------------------------------------#
  #    As linhas de codigo de abaixo sao para poder obter o chute inicial para os graus de liberdade para a ST, SSL e SCN
  #-----------------------------------------------------------------------------------------------------------------------#
  ########################################################################################################################
  
  if(searchNU ==TRUE)
  {
    if(family == "Skew.t")
    {
      mu <- delta <- Delta <- mu1 <- list()
   
      nu               <- seq(lower,upper,space)
      vero             <- c()
      for(i in 1:length(nu))
      {                        
        k1             <- sqrt(nu[i]/2)*gamma((nu[i]-1)/2)/gamma(nu[i]/2)
        b              <- -sqrt(2/pi)*k1
        for (k in 1:g)
        {
          delta[[k]]   <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
          Delta[[k]]   <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
          mu[[k]]      <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
          mu1[[k]]     <- t(matrix(t(x%*%betas[,k]),p,n))
        }
        
        vero[i]    <- sum(log(d.mixedmvST(y, pii, mu, Sigma, shape, nu[i])))## log-likelihood
        if(printNU==TRUE)
          cat(i,"to",length(nu),"\n")
      }
      m           <- cbind(nu,vero)
      maxi        <- max(vero)
      for(i in 1:length(vero))
      {
        if(maxi==m[i,2])
          nu0     <- m[i,1]
      }
      
      if(plotLog == TRUE && saveFigure == TRUE)
      {
        postscript("plotTnu.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
        dev.off() #Fechando o dispositivo potscript
      }
      
      if(plotLog == TRUE  && saveFigure == FALSE)
      {
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
      }
    }
    
    if(family == "Skew.slash")
    {
      mu <- delta <- Delta <- mu1 <- list()
      
      nu               <- seq(lower,upper,space)
      vero             <- c()
      for(i in 1:length(nu))
      {                        
        k1             <- 2*nu[i]/(2*nu[i]-1)
        b              <- -sqrt(2/pi)*k1
        for (k in 1:g)
        {
          delta[[k]]   <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
          Delta[[k]]   <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
          mu[[k]]      <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
          mu1[[k]]     <- t(matrix(t(x%*%betas[,k]),p,n))
        }
        
        vero[i]    <- sum(log(d.mixedmvSS(y, pii, mu, Sigma, shape, nu[i])))## log-likelihood
        if(printNU==TRUE)
          cat(i,"to",length(nu),"\n")
      }
      m           <- cbind(nu,vero)
      maxi        <- max(vero)
      for(i in 1:length(vero))
      {
        if(maxi==m[i,2])
          nu0     <- m[i,1]
      }
      
      if(plotLog == TRUE && saveFigure == TRUE)
      {
        postscript("plotTnu.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
        dev.off() #Fechando o dispositivo potscript
      }
      
      if(plotLog == TRUE  && saveFigure == FALSE)
      {
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
      }
    }
    
    if(family == "Skew.cn")
    {
      mu <- delta <- Delta <- mu1 <- list()
      
      nuu          <- seq(lower,upper,space)
      comp         <- expand.grid(nuu,nuu)
      nu           <- cbind(comp[,1],comp[,2])
      vero         <- c()
      for(i in 1:nrow(nu))
      {                       
        k1             <- nu[i,1]/nu[i,2]^(1/2)+1-nu[i,1]
        b              <- -sqrt(2/pi)*k1
        for (k in 1:g)
        {
          delta[[k]]   <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
          Delta[[k]]   <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
          mu[[k]]      <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
        }
        
        vero[i]    <- sum(log(d.mixedmvSNC(y, pii, mu, Sigma,shape,nu[i,])))
        if(printNU==TRUE)
          cat(i,"to",nrow(nu),"\n")
      }
      m       <- cbind(nu,vero)
      maxi    <- max(vero)
      nu1     <- 0
      nu2     <- 0
      for(i in 1:length(vero))
      {
        if(maxi==m[i,3])
        {
          nu1 <- m[i,1]
          nu2 <- m[i,2]
        }
      }
      nu0     <- c(nu1,nu2)
      
      if(plotLog == TRUE && saveFigure == TRUE)
      {
        z <- matrix(0,nr=length(nuu),nc=length(nuu))
        for(i in 1:length(nuu))
          for(j in 1:length(nuu))
            z[i,j] <- sum(log(d.mixedmvSNC(y, pii, mu, Sigma,shape,c(nuu[i],nuu[j]))))
          postscript("persp.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
          res1 <- persp(x=nuu, y=nuu, z,theta = 95, phi = 15,xlab="nu",ylab="gamma")
          dev.off() #Fechando o dispositivo potscript
          
          postscript("NuGamma.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
          plot(nu[,1],nu[,2],type="n",ylab=expression(paste(gamma)),xlab=expression(paste(nu)),ylim=c(0,1.2))
          abline(v=nu1,lty=2)
          abline(h=nu2,lty=2)
          #eq <- bquote(bold(nu[max] == .(nu0)),bold(nu[max] == .(nu0)))
          #mtext(eq,side=3,cex=1.5)
          legend('topleft',legend= c(as.expression(bquote(nu[max] == .(nu1))), as.expression(bquote(gamma[max] == .(nu2)))),bty="n",cex=1.5)
          dev.off() #Fechando o dispositivo potscript
      }
      
      if(plotLog == TRUE && saveFigure == FALSE)
      {
        z <- matrix(0,nr=length(nuu),nc=length(nuu))
        for(i in 1:length(nuu))
          for(j in 1:length(nuu))
            z[i,j] <- sum(log(d.mixedmvSNC(y, pii, mu, Sigma,shape,c(nuu[i],nuu[j]))))
          par(mfrow=c(1,2))
          res1 <- persp(x=nuu, y=nuu, z,theta = 95, phi = 15,xlab="nu",ylab="gamma")
          plot(nu[,1],nu[,2],type="n",ylab=expression(paste(gamma)),xlab=expression(paste(nu)),ylim=c(0,1.2))
          abline(v=nu1,lty=2)
          abline(h=nu2,lty=2)
          legend('topright',legend= c(as.expression(bquote(nu[max] == .(nu1))), as.expression(bquote(gamma[max] == .(nu2)))),bty="n",cex=1.5)
          par(mfrow=c(1,1))
      }
    }
    
    ###############################################################################################################
    #Objetos de saida
    ###############################################################################################################
    if(family == "Skew.normal"){
      nu0     <- 0
      obj.out <- list(beta = betas, Sigma = Sigma, shape = shape, pii = pii)
    }else{
      obj.out <- list(vero=vero, beta = betas, Sigma = Sigma, shape = shape, pii = pii, nu=nu0)
    }
  }else{
    if(family == "Skew.normal"){
      nu0     <- 0
      obj.out <- list(beta = betas, Sigma = Sigma, shape = shape, pii = pii,nu=nu0)
    }else{
      obj.out <- list(vero=vero, beta = betas, Sigma = Sigma, shape = shape, pii = pii)
    }
    
  }
}

################################################################################
# Function for generate random variables for the class SMSN 
################################################################################
gen.SN.multi <- function(n, mu, Sigma, shape, nu=NULL)
{
 p     <- length(mu)
 delta <- shape / (sqrt(1 + t(shape)%*%shape))
 y     <- matrix(0,n,p)
 for (i in 1:n) 
 {
  y[i,] <- mu + matrix.sqrt(Sigma)%*%(delta*abs(rnorm(1)) + matrix.sqrt(diag(p) - delta%*%t(delta))%*%as.vector(rmvnorm(1, mean = rep(0,p), sigma = diag(p))))
 }
  return(y)
}

gen.ST.multi <- function(n, mu, Sigma, shape, nu)
{
 nu <- as.numeric(nu)
 y  <- matrix(rep(mu, n), n, length(mu), byrow = T) + (rgamma(n, nu/2, nu/2))^(-1/2)*gen.SN.multi(n, rep(0, length(mu)), Sigma, shape)
 return(y)
}

gen.SCN.multi  <- function(n, mu, Sigma, shape, nu)
{
 x1 <- matrix(0,n,length(mu))
 for (i in 1:n)
 {
  u <- runif(1)
  if (u < nu[1]) {x1[i,] <- gen.SN.multi(1, mu, Sigma/nu[1], shape)}
  if (u > nu[1]) {x1[i,] <- gen.SN.multi(1, mu, Sigma, shape)}
 }
 return(x1)
}

gen.SS.multi  <- function(n, mu, Sigma, shape, nu)
{
 u1 <- runif(n)
 u2 <- u1^(1/(nu))   # formula 10 do artigo e metodo da inversao
 ys <- mu + (u2)^(-1/2)*gen.SN.multi(n, c(0,0), Sigma, shape)
 return(ys)
}
################################################################################

dadosSimul <- function(n, pii, beta, Sigma, shape, nu, family, a=-1, b=1)
{
 if(sum(pii) != 1) stop("Sum of pii diferent of one")
 if(a >= b) stop("Limits of runif of X, a > b") 
 
 w  <- rmultinom(n, size = 1, prob = pii)
 G  <- rowSums(w)
 z  <- h <- w <- NULL
 q  <- length(beta[[1]])
 p  <- ncol(Sigma[[1]])
  
 mu <- list()
  
 for (r in 1:length(G))
 {
   # mus[[r]] <- X%*%beta[[r]]
  X        <- NULL
  for(s in 1:(q-1)) X <- cbind(X,runif(G[r]*p,a,b))
  X        <- cbind(1,X)
  mu[[r]]  <- matrix(X%*%beta[[r]], ncol = p, byrow = TRUE)
  y        <- matrix(0, G[r], p)
    
  for(l in 1:G[r])
  {
   if(family[[r]] == "Normal")
   {
    shape[[r]] <- rep(0,p)
    y[l,]      <- gen.SN.multi(n = 1, mu = mu[[r]][l,], Sigma = Sigma[[r]], shape = shape[[r]]) 
   }
      
   if(family[[r]] == "t")
   {
    shape[[r]] <- rep(0,p)
    y[l,]      <- gen.ST.multi(n = 1, mu = mu[[r]][l,], Sigma = Sigma[[r]], shape = shape[[r]], nu = nu[r])
   } 
   if(family[[r]] == "Skew.normal") y[l,] <- gen.SN.multi(n = 1, mu = mu[[r]][l,], Sigma = Sigma[[r]], shape = shape[[r]], nu = nu[r])
   if(family[[r]] == "Skew.t")      y[l,] <- gen.ST.multi(n = 1, mu = mu[[r]][l,], Sigma = Sigma[[r]], shape = shape[[r]], nu = nu[r])
   if(family[[r]] == "Skew.cn")     y[l,] <- gen.SCN.multi(n = 1, mu = mu[[r]][l,], Sigma = Sigma[[r]], shape = shape[[r]], nu = nu[r])
   if(family[[r]] == "Skew.slash")  y[l,] <- gen.SS.multi(n = 1, mu = mu[[r]][l,], Sigma = Sigma[[r]], shape = shape[[r]], nu = nu[r])   
  } # fim for l
    
    m   <- dim(y)[1]
    y <- t(t(y))
    z <- rbind(z,y)
    h <- rbind(h,cc)
    w <- rbind(w,X)
    
  } # fim for G  
  
  return(list(y = z, X = w, cc = h, G = G))
}


d.mixed.multi.SMSN.Reg <- function(y, pi1, xB, Sigma){
  #y: Matrix de dados m x p
  #pi1: deve ser do tipo vetor de dimensao g
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
  
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dmvNCensReg(cc, y , xB[[j]], Sigma[[j]])
  return(dens)
}
 
#################################################################################################################################
#model=TrySMSN
#model=fitSN_g2
imm.smsn.mmrm <- function(y, x, model)
{
  adjoint     <- function(A) det(A)*solve(A)
  deriv.der   <- function(A,B,C) det(A)*sum(B * t(C))
  
  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
  if (ncol(y) <= 1) stop(paste("The dimension of y (p) is: ", ncol(y),". We need p >= 2.",sep=" "))
  #if (model$uni.Gama) stop("Sorry. The information matrix cannot be calculated when the uni.Gama was used!")
  y           <- as.matrix(y)
  dimnames(y) <- NULL
  n           <- nrow(y)
  p           <- ncol(y)
  g           <- length(model$pii)
  Sipi        <- Simu <- Silambda <- c()
  Ssigma      <- c()
  Sbeta       <- c()
  
  for(i in 1:length(model$Sigma)) model$Sigma[[i]] <- model$Sigma[[i]] %*% model$Sigma[[i]]    
  
  betas       <- model$beta
  
  q           <- ncol(x)
  
  Sigma       <- model$Sigma
  shape       <- model$shape
  pii         <- model$pii
  nu          <- model$nu
  
  if (class(model) == "Skew.t")
  {
    k1               <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    b                <- -sqrt(2/pi)*k1
    
    Delta <- delta <- mu <-  list()
    
    for(k in 1:g)
    {
      delta[[k]]     <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
      Delta[[k]]     <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
      mu[[k]]        <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
    }
    
    soma             <- soma2 <- 0
    
    I.Phi            <- function(w=0,Ai=NULL,di,nu=0) as.numeric((( 2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt( ((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))
    I.phi            <- function(w=0,Ai=NULL,di,nu=0) as.numeric(((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))
    nj               <- rep(p,n)
    
    for (i in 1:n)
    {
      S              <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu       <- 0
      yi             <- matrix(y[i,], 1, p)
      
      for (j in 1:g)
      {
        
        x1           <- matrix(x[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ], ncol=q)  
        Dr           <- matrix.sqrt(Sigma[[j]])
        Dr.inv       <- solve(Dr)
        d.sig        <- det(Sigma[[j]])        
        
        Ai           <- as.numeric(t(shape[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]][i,]))
        di           <- as.numeric(mahalanobis(yi, mu[[j]][i,], Sigma[[j]]))
        
        dir.dbeta    <- -2*t(x1)%*%(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]][i,])
        dAir.dbeta   <- -t(x1)%*%Dr.inv%*%shape[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]][i,])
        
        dPsi.dbeta   <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dbeta * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dbeta*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*  dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)
        
        #para os elementos de sigma                      
        l <- m <- 1
        
        for(k in 1:((p+1)*p/2))
        {
          Vis         <- FALSE
          D           <- matrix(rep(0,p*p),p,p)
          D[l,m]      <- D[m,l] <- 1
          
          ddet.ds     <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
          dir.ds      <- - t(y[i,] - mu[[j]][i,])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          dAir.ds     <- - t(shape[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          
          dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
          Ssigma[k]   <- (pii[j]/ d.mixedmvST(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu) )*dPsi.dsigma
          
          
          if(((l*m - p*floor((l*m)/p)) == 0) && (l != m))
          {
            l    <- l+1
            m    <- l
            Vis  <- TRUE
          }
          
          if(!Vis) m <- m+1
        }
        
        ui       <- rgamma(10000, shape = nu/2, rate = nu/2)
        resto    <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        dPsi.dnu <- dPsi.dnu + pii[j]*((d.sig^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)
        
        Simu     <- as.vector((pii[j]/ d.mixedmvST(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu) )*dPsi.dbeta)
        Silambda <- as.vector((pii[j]/ d.mixedmvST(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu) )*dPsi.dlambda )
        
        S        <- c(S, Simu, Silambda, Ssigma)
      }
      
      Sinu       <- (1/d.mixedmvST(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu))*dPsi.dnu
      
      if(g>1)
      {
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvST(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu))*( dmvt.ls(yi, mu[[j]][i,], Sigma[[j]], shape[[j]], nu) - dmvt.ls(yi, mu[[g]][i,], Sigma[[g]], shape[[g]], nu))        
        S        <- c(S, Sipi, Sinu)
      }
      
      if(g == 1) {S <- c(S, Sinu)}
      soma       <- soma + S%*%t(S)
    }

    
  }
  
  if (class(model) == "Skew.cn")
  {
    k1          <- nu[1]/nu[2]^(1/2)+1-nu[1]
    b           <- -sqrt(2/pi)*k1
    
    Delta <- delta <- mu <-  list()
    
    for(k in 1:g)
    {
      delta[[k]]     <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
      Delta[[k]]     <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
      #mu[[k]]        <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
      mu[[k]]        <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
    }
    
    soma        <- soma2 <- 0
    
    I.Phi       <- function(w=0,Ai=NULL,di,nu=0) as.numeric( sqrt(2*pi)*(nu[1]*nu[2]^(w -0.5)*dnorm(sqrt(di), 0, sqrt(1/nu[2]))*pnorm(nu[2]^(1/2)*Ai) + (1 - nu[1])*(dnorm(sqrt(di), 0,1)*pnorm(Ai)) )   )
    I.phi       <- function(w=0,Ai=NULL,di,nu=0) as.numeric( nu[1]*nu[2]^(w - 0.5)*dnorm(sqrt(di + Ai^2), 0, sqrt(1/nu[2])) + (1 - nu[1])*dnorm(sqrt(di + Ai^2))   )
    nj          <- rep(p,n)
    
    for (i in 1:n)
    {
      S         <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu1 <- dPsi.dnu2 <- 0
      yi        <- matrix(y[i,], 1, p)#2)
      x1        <- matrix(x[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ], ncol=q)  
      
      for (j in 1:g)
      {
        Dr           <- matrix.sqrt(Sigma[[j]])
        Dr.inv       <- solve(Dr)
        d.sig        <- det(Sigma[[j]])        
        
        Ai           <- as.numeric(t(shape[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]][i,]))
        di           <- as.numeric(mahalanobis(yi, mu[[j]][i,], Sigma[[j]]))
        
        #derivadinhas
        dir.dbeta    <- -2*t(x1)%*%(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]][i,])
        dAir.dbeta   <- -t(x1)%*%Dr.inv%*%shape[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]][i,])
        
        dPsi.dbeta   <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dbeta * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dbeta*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)
        
        dPsi.dnu1    <- dPsi.dnu1 + pii[j]*2*(dmvnorm(yi, mu[[j]][i,], nu[2]^(-1)*Sigma[[j]])*pnorm(nu[2]^(1/2)*Ai) - dmvnorm(yi, mu[[j]][i,], Sigma[[j]])*pnorm(Ai) )
        dPsi.dnu2    <- dPsi.dnu1 + pii[j]*((nu[1]*d.sig^(-1/2)*nu[2]^(p/2))/(2*pi)^(p/2))*exp(-nu[2]*di/2)*(p*nu[2]^(-1)*pnorm(nu[2]^(1/2)*Ai) + dnorm(nu[2]^(1/2)*Ai)*Ai*nu[2]^(-1/2) - pnorm(nu[2]^(1/2)*Ai)*di )
        
        #para os elementos de sigma
        l   <- m <- 1
        for(k in 1:((p+1)*p/2))
        {
          Vis         <- FALSE
          D           <- matrix(rep(0,p*p),p,p)
          D[l,m]      <- D[m,l] <- 1
          
          ddet.ds     <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
          dir.ds      <- - t(y[i,] - mu[[j]][i,])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          dAir.ds     <- - t(shape[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          
          dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
          Ssigma[k]   <- (pii[j]/ d.mixedmvSNC(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu) )*dPsi.dsigma
          
          if(((l*m - p*floor((l*m)/p)) == 0) && (l != m))
          {
            l   <- l+1
            m   <- l
            Vis <- TRUE
          }
          if(!Vis) m <- m+1
        }
        
        
        Simu     <- as.vector((pii[j]/ d.mixedmvSNC(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu) )*dPsi.dbeta)
        Silambda <- as.vector((pii[j]/ d.mixedmvSNC(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu) )*dPsi.dlambda )
        
        S        <- c(S, Simu, Silambda, Ssigma)
      }
      Sinu1 <- (1/d.mixedmvSNC(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu))*dPsi.dnu1
      Sinu2 <- (1/d.mixedmvSNC(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu))*dPsi.dnu2
      
      if (g >1)
      { 
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSNC(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape, nu))*( dmvSNC(yi, mu[[j]][i,], Sigma[[j]], shape[[j]], nu) - dmvSNC(yi, mu[[g]][i,], Sigma[[g]], shape[[g]], nu))
        S <- c(S, Sipi, Sinu1, Sinu2)
      }
      if( g==1) S <- c(S, Sinu1, Sinu2)
      soma <- soma + S%*%t(S)
      #      soma2 <- soma2 + S
    }
    
#    NAME <- piiN <- c()
#    for(i in 1:g){
#      SigmaN <- muN <- shapeN <- c()
    #      for (k in 1:p){
    #        muN <- c(muN,paste("mu",i,"_",k,sep=""))
    #        shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
    #      }
    #    l <- m <- 1
    #  for(k in 1:((p+1)*p/2)) {
    #    Vis <- FALSE
    #    SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
    #    if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
    #      l <- l+1
    #      m <- l
    #      Vis <- TRUE
    #    }
    #    if(!Vis) m <- m+1
    #  }
    #  NAME <- c(NAME,muN,shapeN,SigmaN)
  # }
    
  #  if( g==1) NAME <- c(NAME,"nu1","nu2")
  #  else{
  #    for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
  #    NAME <- c(NAME,piiN,"nu1","nu2")
  #  }
    
   #  dimnames(soma)[[1]] <- NAME
   # dimnames(soma)[[2]] <- NAME   
    
  }
  
  if (class(model) == "Skew.slash")
  {
    k1  <- 2*nu/(2*nu-1)
    b   <- -sqrt(2/pi)*k1
    
    Delta <- delta <- mu <-  list()
    
    for(k in 1:g)
    {
      delta[[k]]     <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
      Delta[[k]]     <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
      mu[[k]]        <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
    }
    
    soma <- soma2 <- 0
    
    I.Phi <- function(w=0,Ai=NULL,di,nu=0)
    {
      Esper <- vector(mode = "numeric", length = length(di))
      for(i in 1:length(di))
      {
        U        <- runif(2500)
        V        <- pgamma(1,w + nu, di[i]/2)*U
        S        <- qgamma(V,w + nu, di[i]/2)
        Esper[i] <- mean(pnorm(S^(1/2)*Ai[i]))
      }
      res1       <- (nu*(2^(w + nu)*gamma(w + nu))/(di^(w + nu)))*pgamma(1, w + nu, di/2)*Esper
      return(res1)
    }
    
    I.phi <- function(w=0,Ai=NULL,di,nu=0)
    {
     res2 <- ((nu*2^(w + nu)*gamma(w + nu))/(sqrt(2*pi)*(di + Ai^2)^(w + nu)))*pgamma(1, w + nu, (di + Ai^2)/2)
     return(res2)
    }
    
    nj         <- rep(p,n)
    
    
    for (i in 1:n)
    {
      S        <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      yi       <- matrix(y[i,], 1, p)#2)
      x1       <- matrix(x[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ], ncol=q)  
      for (j in 1:g)
      {
        Dr     <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig  <- det(Sigma[[j]])        
        
        Ai     <- as.numeric(t(shape[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]][i,]))
        di     <- as.numeric(mahalanobis(yi, mu[[j]][i,], Sigma[[j]]))
        
        #derivadinhas
        dir.dbeta    <- -2*t(x1)%*%(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]][i,])
        dAir.dbeta   <- -t(x1)%*%Dr.inv%*%shape[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]][i,])
        
        dPsi.dbeta   <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dbeta * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dbeta*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)
        
        #para os elementos de sigma
        l <- m <- 1
        for(k in 1:((p+1)*p/2))
        {
          Vis         <- FALSE
          D           <- matrix(rep(0,p*p),p,p)
          D[l,m]      <- D[m,l] <- 1
          
          ddet.ds     <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
          dir.ds      <- - t(y[i,] - mu[[j]][i,])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          dAir.ds     <- - t(shape[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          mu_2        <- do.call("list", lapply(mu, "[", i, ))
          mu_2[[1]]   <- rbind( mu_2[[1]])
          mu_2[[2]]   <- rbind( mu_2[[2]])
          
          dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
          Ssigma[k]   <- (pii[j]/ d.mixedmvSS(yi, pii, mu_2, Sigma, shape, nu) )*dPsi.dsigma
          
          if(((l*m - p*floor((l*m)/p)) == 0) && (l != m))
          {
            l   <- l+1
            m   <- l
            Vis <- TRUE
          }
          if(!Vis) m <- m+1
        }
        
        #ui <- rgamma(10000, shape = nu/2, rate = nu/2)
        #resto <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        #dPsi.dnu <- dPsi.dnu + pii[j]*((det(Sigma[[j]])^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-4*digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)
        
        u        <- runif(8000)
        dPsi.dnu <- dPsi.dnu + pii[j]*mean(2*u^(nu - 1)*(1 + nu*log(u))*( (det(Sigma[[j]])^(-1/2)/(2*pi)^(p/2))*exp(-(u^(-1)/2)*di) )*pnorm(u^(1/2)*as.numeric(t(shape[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%t(yi-mu[[j]][i,]))))#integrate(f,0,1)$value
        
        Simu     <- as.vector((pii[j]/ d.mixedmvSS(yi, pii, mu_2, Sigma, shape, nu) )*dPsi.dbeta)
        Silambda <- as.vector((pii[j]/ d.mixedmvSS(yi, pii, mu_2, Sigma, shape, nu) )*dPsi.dlambda )
        
        S        <- c(S, Simu, Silambda, Ssigma)
      }
      Sinu <- (1/d.mixedmvSS(yi, pii, mu_2, Sigma, shape, nu))*dPsi.dnu
      if(g>1){
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSS(yi, pii, mu_2, Sigma, shape, nu))*(dmvSS(yi,  rbind(mu[[j]][i,]), Sigma[[j]], shape[[j]], nu) - dmvSS(yi, rbind(mu[[g]][i,]), Sigma[[g]], shape[[g]], nu))
        S <- c(S, Sipi, Sinu)                                          
      }
      if( g==1) S <- c(S, Sinu)
      soma <- soma + S%*%t(S)
      #      soma2 <- soma2 + S
    }
#    NAME <- piiN <- c()
#    for(i in 1:g){
#      SigmaN <- muN <- shapeN <- c()
#      for (k in 1:p){
#        muN <- c(muN,paste("mu",i,"_",k,sep=""))
#        shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
#      }
#      l <- m <- 1
#      for(k in 1:((p+1)*p/2)) {
#        Vis <- FALSE
#        SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
#        if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
#          l <- l+1
#          m <- l
#          Vis <- TRUE
#        }
#        if(!Vis) m <- m+1
#      }
#      NAME <- c(NAME,muN,shapeN,SigmaN)
#    }
    
#    if( g==1) NAME <- c(NAME,"nu")    
#    else{
#      for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
#      NAME <- c(NAME,piiN,"nu")
#    }
#    
#    dimnames(soma)[[1]] <- NAME
#    dimnames(soma)[[2]] <- NAME   
  }
  
  if (class(model) == "Skew.normal")
  {
    k1               <- 1
    b                <- -sqrt(2/pi)*k1
    
    Delta <- delta <- mu <-  list()
    
    for(k in 1:g)
    {
      delta[[k]]     <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
      Delta[[k]]     <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
      mu[[k]]        <- t(matrix(t(x%*%betas[,k]),p,n)) + matrix(rep(b*Delta[[k]], n), n, p, byrow = TRUE)
    }
    
    soma             <- soma2 <- 0
    
    I.Phi            <- function(w=0,Ai=NULL,di) as.numeric(exp(-di/2)*pnorm(Ai))   
    I.phi            <- function(w=0,Ai=NULL,di)   as.numeric(exp(-di/2)*dnorm(Ai))       
    nj               <- rep(p,n)
    
    for (i in 1:n)
    {
      S              <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu       <- 0
      yi             <- matrix(y[i,], 1, p)
      x1             <- matrix(x[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),  ], ncol=q)  
      
      for (j in 1:g)
      {
        Dr           <- matrix.sqrt(Sigma[[j]])
        Dr.inv       <- solve(Dr)
        d.sig        <- det(Sigma[[j]])        
        
        Ai           <- as.numeric(t(shape[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]][i,]))
        di           <- as.numeric(mahalanobis(yi, mu[[j]][i,], Sigma[[j]]))
        
        dir.dbeta    <- -2*t(x1)%*%(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]][i,])
        dAir.dbeta   <- -t(x1)%*%Dr.inv%*%shape[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]][i,])
        
        dPsi.dbeta   <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dbeta * I.phi((p+1)/2, Ai, di) - (1/2)*dir.dbeta*I.Phi((p/2)+1, Ai, di) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*  dAir.dlambda*I.phi((p+1)/2, Ai, di)
        
        #para os elementos de sigma                      
        l <- m <- 1
        
        for(k in 1:((p+1)*p/2))
        {
          Vis         <- FALSE
          D           <- matrix(rep(0,p*p),p,p)
          D[l,m]      <- D[m,l] <- 1
          
          ddet.ds     <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
          dir.ds      <- - t(y[i,] - mu[[j]][i,])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          dAir.ds     <- - t(shape[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]][i,])
          
          dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di) )           
          Ssigma[k]   <- (pii[j]/d.mixedmvSN(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape))*dPsi.dsigma
          
          
          if(((l*m - p*floor((l*m)/p)) == 0) && (l != m))
          {
            l    <- l+1
            m    <- l
            Vis  <- TRUE
          }
          
          if(!Vis) m <- m+1
        }
        
        Simu     <- as.vector((pii[j]/ d.mixedmvSN(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape) )*dPsi.dbeta)
        Silambda <- as.vector((pii[j]/ d.mixedmvSN(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape) )*dPsi.dlambda )
        
        S        <- c(S, Simu, Silambda, Ssigma)
      }
      
      if(g>1)
      {
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSN(yi, pii, do.call("list", lapply(mu, "[", i, )), Sigma, shape))*( dmvSN(yi, mu[[j]][i,], Sigma[[j]], shape[[j]]) - dmvSN(yi, mu[[g]][i,], Sigma[[g]], shape[[g]]))        
        S        <- c(S, Sipi)
      }
      
      if(g == 1) {S <- c(S)}
      soma       <- soma + S%*%t(S)
    }
  }
  
  return(list(IM=soma, EP=sqrt(diag(solve(soma)))))  
  #  return(list(soma, soma2))
}
#EP <- imm.smsn.mmrm(y, x, fitSN_g2)
#EP$EP


