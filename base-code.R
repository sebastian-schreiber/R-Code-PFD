# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This file provides the  base commands 
# for simulating and exploring the deterministic, two species competition model 
# Ni[t+1]=Ni[t]*(Ni[t]/(ai+Ni[t]+bi*Nj[t]))*(lambdai/(1+alphai*Ni[t]+betai*Nj[t]))
# where lambdai is the intrinsic growth rate of species i
# alphai/betai are the strength of intra/interspecific competition on species i
# ai and bi determine the strength of positive frequency-dependence
# Most of this code (as with the Ecology paper) only examine the model with ai=0

# required packages
require(rootSolve)



# The following code creates figures with 
# nullclines and sample solutions
# input: list of parms lambda, a, b, alpha, beta
# output: plot of nullclines and sample solutions 
nullclines=function(parms=parms,main.text=""){
    with(as.list(parms),{
      f1=function(x,y)lambda[1]*x-(a[1]+x+b[1]*y)*(1+alpha[1]*x+beta[1]*y)
      f2=function(x,y)lambda[2]*y-(a[2]+y+b[2]*x)*(1+alpha[2]*y+beta[2]*x)
      L=50
      xs=seq(0,M,length=L)
      ys=seq(0,M,length=L)
      M1=outer(X = xs,Y = ys,FUN = f1)
      contour(xs,ys,M1,levels=c(0,0),col="red",lwd=2,labcex=1,drawlabels=FALSE,xlab="species 1 density",ylab="species 2 density",main=main.text)
      M2=outer(X = xs,Y = ys,FUN = f2)
      contour(xs,ys,M2,levels=c(0,0),add=TRUE,col="blue",lwd=2,labcex=1,drawlabels =FALSE)
      #legend("topright",c("species 1 nullcline","species 2 nullcline"),col=c("red","blue"),lty=1,lwd=2,bty="n") # deterministic version
      legend("topright",c("species 1 nullclines","species 2 nullclines","community trajectory"),col=c("red","blue","gray"),lty=1,lwd=2,bty="n") # stochastic verison
      Tf=100
      X=matrix(NA,Tf,2)
      G=function(x)c(x[1]*lambda[1]*x[1]/(a[1]+x[1]+b[1]*x[2])/(1+alpha[1]*x[1]+beta[1]*x[2]),
                     x[2]*lambda[2]*x[2]/(a[2]+x[2]+b[2]*x[1])/(1+alpha[2]*x[2]+beta[2]*x[1]) )
      if(length(X0)>0){for(j in 1:length(X0[,1])){
        X[1,]=X0[j,]
        for(t in 2:Tf)X[t,]=G(X[t-1,])
        lines(X[,1],X[,2],col="gray",cex=0.5)
        if(is.element(j,Xarrows)){
          for(t in 2:5){
          dx=X[t+1,1]-X[t,1]
          dy=X[t+1,2]-X[t,2]
          arrows(x0=X[t,1],y0=X[t,2],x1=X[t,1]+dx/4,y1=X[t,2]+dy/4,length = 0.1,angle = 20 ,col="gray")
          }}
          }}
      contour(xs,ys,M1,levels=c(0,0),col="red",lwd=2,labcex=1,drawlabels=FALSE,add=TRUE)
      contour(xs,ys,M2,levels=c(0,0),col="blue",lwd=2,labcex=1,drawlabels =FALSE,add=TRUE)
    }
    )
  }
  
# The following code solves for 
# equilibria given an initial guess
# input: list of parms lambda, a, b, alpha, beta and initial guess
# output: plot of the equilibrium of the specified color

equi=function(parms,guess,coll="black",radial=FALSE){
    with(as.list(parms),{
      f1=function(x,y)lambda[1]*x^2/(a[1]+x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
      f2=function(x,y)lambda[2]*y^2/(a[2]+y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
      h=function(x){(f1(x[1],x[2])-x[1])^2+(f2(x[1],x[2])-x[2])^2}
      out=optim(guess,h)
      if(radial)lines(x=c(0,out$par[1])*3,y=c(0,out$par[2])*3,lty=2)
      points(x=out$par[1],y=out$par[2],pch=21,bg=coll,cex=1.75)
    }
    )
  }
  
# the following code   
# numerically computes the carrying simplex
# input: list of parms lambda, a, b, alpha, beta
# output: vectors of x and y values corresponding to the carrying simplex
 
  carrying=function(parms){
    with(as.list(parms),{
      #assuming that a[i]=0
      f1=function(x,y)lambda[1]*x^2/(a[1]+x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
      f2=function(x,y)lambda[2]*y^2/(a[2]+y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
      L=10000
      eq<-(lambda-1)/alpha
      c<-seq(0,1,length=L)
      xs<-c*eq[1]
      ys<-(1-c)*eq[2]
      Tf<-50
      for(t in 1:Tf){
        xs.new<-f1(xs,ys)
        ys.new<-f2(xs,ys)
        xs<-xs.new
        ys<-ys.new
      }
      return(cbind(xs,ys))
    }
    )
  }
  
  
# the following code 
# computes relative fitness across equilibrium transect 
# input: list of parms lambda, a, b, alpha, beta
# output plots of the relative fitnesses
  
  fit=function(parms){
    with(as.list(parms),{
      f1=function(x,y)lambda[1]*x/(a[1]+x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
      f2=function(x,y)lambda[2]*y/(a[2]+y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
      xstar<-(lambda[1]-1)/alpha[1]
      ystar<-(lambda[2]-1)/alpha[2]
      g1=function(c)f1(c*xstar,(1-c)*xstar)
      g2=function(c)f2(c*xstar,(1-c)*xstar)
      cs=seq(0,1,length=100)
      f1s=g1(cs)
      f2s=g2(cs)
      matplot(cs,cbind(f1s/f2s,f2s/f1s),type="l",xlab="species 1 frequency",
              ylab=expression(paste("relative fitness ",f[1]/f[2])),lwd=3,lty=1,bty="n",
              col=c("red","blue"),ylim=c(0,2))
      abline(h=1,lty=2,lwd=3)
      h1=function(c)g1(c)-g2(c)
      out=uniroot.all(h1,c(0,1))
      temp=which(out*(1-out)>0)
      abline(v=max(out[temp]),lty=2)
      abline(v=min(out[temp]),lty=2)
    }
    )
  }
  
  
  
  
# the following code 
# computes relative fitness across equilibrium transect
# input: list of parameters, leg whether or not to include a legend in plot, approx whether or not to show the approximation
# output plots of the relative fitnesses  

    logfit=function(parms,leg=TRUE,approx=TRUE){
    with(as.list(parms),{
      f1=function(x,y)lambda[1]*x/(a[1]+x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
      f2=function(x,y)lambda[2]*y/(a[2]+y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
      RI1=function(x,y)x/(a[1]+x+b[1]*y)
      RI2=function(x,y)y/(a[2]+y+b[2]*x)
      xstar<-(lambda[1]-1)/alpha[1]
      ystar<-(lambda[2]-1)/alpha[2]
      g1=function(c)f1(c*xstar,(1-c)*xstar)
      g2=function(c)f2(c*xstar,(1-c)*xstar)
      h1=function(c)RI1(c*xstar,(1-c)*xstar)
      h2=function(c)RI2(c*xstar,(1-c)*xstar)
      cs=seq(0.0001,0.9999,length=100)
      f1s=g1(cs)
      f2s=g2(cs)
      RI1s=h1(cs)
      RI2s=h2(cs)
      tempcols=c("black","black","blue")
      templtys=c(2,3,1)
      matplot(cs,cbind((f1s/RI1s)/(f2s/RI2s),(RI1s)/(RI2s),f1s/f2s),type="l",xlab="species 1 frequency",
              ylab=expression(paste("")),lwd=3,bty="n",
              col=tempcols,log="y",lty=templtys,ylim=c(0.05,5))
      abline(h=1,lty=1,lwd=1,col="darkgray")
      if(approx){
        f1=function(x,y)lambda[1]*x/(a[1]+x+b[1]*y)/(alpha[1]*x+beta[1]*y)
        f2=function(x,y)lambda[2]*y/(a[2]+y+b[2]*x)/(alpha[2]*y+beta[2]*x)
        g1=function(c)f1(c*xstar,(1-c)*xstar)
        g2=function(c)f2(c*xstar,(1-c)*xstar)
        f1s=g1(cs)
        f2s=g2(cs)
        matplot(cs,cbind((f1s/RI1s)/(f2s/RI2s),f1s/f2s),type="l",col=rgb(1,0,0,0.3),lwd=6,lty=1,add=TRUE)
      }
      h1=function(c)g1(c)-g2(c)
      out=uniroot.all(h1,c(0,1))
      temp=which(out*(1-out)>0)
      eqs=sort(out[temp])
      cols=c("white","black","white")
      for(j in 1:length(temp)){
        points(eqs[j],1,pch=21,bg=cols[j],cex=1.75)
      }
      if(leg)legend("bottom",c("relative per-capita growth\n w/o PFD","relative PFD","relative per-capita growth"),lty=templtys,col=tempcols,bty="n",cex=1)
    }
    )
  }
  
  

  
# the following code  
# computes the roots for the highly productive limit (i.e. lambdai XXL)
# the calculations for the relevant polynomial are in 
# the supplement S1 
# input: list of parms lambda, a, b, alpha, beta
# output: equilibria for the ratio dynamics (see supplement S1)
  
high=function(parms){
    with(as.list(parms),{
      z<-numeric(4);
      z[1]<--lambda[2]*b[2]*beta[2]
      z[2]<-lambda[1]*alpha[1]-lambda[2]*(alpha[1]*b[2]+beta[2])
      z[3]<-lambda[1]*(b[1]*alpha[2]+beta[1])-lambda[2]*alpha[1]
      z[4]<-lambda[1]*b[1]*beta[1]
      return(polyroot(z))
    }
    )
  }
  
  
  
# an alternative version of the previous
# command which only returns the real parts 
# to find the roots in the limiting case of 
# highly productive systems
  equi2=function(parms){
    with(as.list(parms),{
      z=numeric(4)
      z[1]=-lambda[2]*beta[1]*b[1]
      z[2]=-lambda[2]*b[1]-lambda[2]*beta[1]+lambda[1]
      z[3]=-lambda[2]+lambda[1]*b[2]+lambda[1]*beta[2]
      z[4]=lambda[1]*b[2]*beta[2]
      zstar=polyroot(z)
      relevant=which((abs(Im(zstar))<10^(-7))&(Re(zstar)>0))
      return(Re(zstar[relevant]))
    }
    )
  }