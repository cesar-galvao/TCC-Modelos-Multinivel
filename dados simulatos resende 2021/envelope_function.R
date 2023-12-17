#*************************************************************************Description*************************************************************************************************#
#Create normal probability plots of residuals with simulated envelope to assess the goodness of normal regression fit under OLS and WOLS estimators
#***ARGUMENTS***#
# fit - an object containing the results returned by lm() fitting.
# OLS - logical; if TRUE, the OLS is used to create the plot; if FALSE, the WOLS is used to create the plot. Default is TRUE.
# main.title - main title for the plot. Default is "Envelope".
# faixa.fixed - range of residuals values (optional). Default is NULL.
# labels.fixed - labels of the observations used to create the plot (optional). Default is NULL meaning that all observations are used.
#**************************************************************************************************************************************************************************************#

envelope_LR <- function(fit, OLS=T, main.title = "Envelope", faixa.fixed = NULL, labels.fixed = NULL) { 
  B <- 100; #number of replicates
  X <- model.matrix(fit)
  p <- ncol(X);  n <- nrow(X)
  
  #***parameters for parametric bootstrap***#
  
  Menvelope_r <- matrix(numeric(0),nrow=n,ncol=B)
  
  #------------> residuals for the observed sample<--------------#
  ts <- rstudent(fit) #Studentized residuals
  betahat <- as.numeric(fit$coefficients)
  sigma2hatc <- sum(resid(fit)^2)/(n-p) #sigma2 estimate (constant)
  
  if(OLS==T){
    sigma2hat <- rep(sigma2hatc,n)
  }
  else sigma2hat <- sigma2hatc/fit$weights
  
  for(j in 1:B){		
    ygen <- rnorm(n, X%*%betahat, sqrt(sigma2hat))
    
    if(OLS==T){
      fitb <- lm(ygen~X[,2:p])
    }
    else{
      fitb <- lm(ygen~X[,2:p], weights= fit$weights)
    }
    
    Res <- rstudent(fitb)
    Menvelope_r[,j] = Res
  }
  Menvelope <- apply(Menvelope_r,2,sort);          
  res <-    ts;    
  res_min  <-    as.numeric(t(apply(t(Menvelope), 2,quantile, probs =0.05)));         
  res_mean <-    as.numeric(t(apply(t(Menvelope), 2,quantile, probs =0.5)));                              
  res_max  <-    as.numeric(t(apply(t(Menvelope), 2,quantile, probs =0.95)));           
  faixa <- range(res,res_min,res_max)
  if(is.vector(faixa.fixed)) faixa <- faixa.fixed
  if(is.vector(labels.fixed)) labels <- labels.fixed
  par(mar=c(5.0,5.0,4.0,2.0))
  v <- qqnorm(res, main=main.title, xlab="Quantis da Normal", ylab="Resíduos Studentizados", ylim=faixa, pch=16, cex=.7, cex.lab=1.2, cex.axis=1.5, cex.main=1.5)
  #identify(v$x,v$y,labels,cex =1.3) #identify points in the plot
  par(new=T)
  #
  qqnorm(res_min,axes=F,main = "",xlab="",ylab="",type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_max,axes=F,main = "",xlab="",ylab="", type="l",ylim=faixa,lty=1,lwd=2.0)
  par(new=T)
  qqnorm(res_mean,axes=F,xlab="",main = "", ylab="", type="l",ylim=faixa,lty=2,lwd=2.0)
}#ends function
