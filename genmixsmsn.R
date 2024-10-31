genmixsmsn = function(n, family, mu1,mu2,Sigma1,Sigma2,shape1,shape2,nu=NULL)
{
 y <- matrix(0,n,2)  
 if(family=="Skew.normal")
  for (i in 1:n){
   arg1   = list(mu=c(mu1[i,]), Sigma=Sigma1, shape=c(shape1)); arg2   = list(mu=c(mu2[i,]), Sigma=Sigma2, shape=c(shape2))
   y[i,]  <- rmmix(1, pii, family, list(arg1,arg2))}
  
 if(family=="Skew.t" || family=="Skew.slash" || family=="Skew.cn")
  for (i in 1:n){
   arg1   = list(mu=c(mu1[i,]), Sigma=Sigma1, shape=c(shape1), nu=nu); arg2   = list(mu=c(mu2[i,]), Sigma=Sigma2, shape=c(shape2), nu=nu)
   y[i,]  <- rmmix(1, pii, family, list(arg1,arg2))}
 
 return(y)
}