source("base-code.R")
cs<-seq(0,1,length=30)
parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=150,
  X0=c())

pdf("why-rho.pdf",width=8,height=4)
#dev.new(width=8,height=4); 
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1), cex.lab = 1.5, cex.axis = 1.5)
nullclines(parms)

source("base-code-stochastic.R")
set.seed(3)
parms=list(b=c(1,1)/4,alpha=c(1,1),beta=c(1,1)/5,
           lambda=c(100,100),sd.log.lambda=c(0.05,0.05),
           rho=-1,reps=2,Tf=55,a=c(0,0))
out=equi(parms,c(100,100))
run=um(parms=parms,x0=out$par)

nullclines1=function(parms){
    with(as.list(parms),{
      f1=function(x,y)lambda[1]*x-(a[1]+x+b[1]*y)*(1+alpha[1]*x+beta[1]*y)
      f2=function(x,y)lambda[2]*y-(a[2]+y+b[2]*x)*(1+alpha[2]*y+beta[2]*x)
      L=50
      xs=seq(0,M,length=L)
      ys=seq(0,M,length=L)
      M1=outer(X = xs,Y = ys,FUN = f1)
      M2=outer(X = xs,Y = ys,FUN = f2)
      contour(xs,ys,M1,levels=c(0,0),col=rgb(1, 0, 0, alpha=0.1),lwd=0.1,labcex=1,drawlabels=FALSE,add=TRUE)
      contour(xs,ys,M2,levels=c(0,0),col=rgb(0, 0, 1, alpha=0.1),lwd=0.1,labcex=1,drawlabels =FALSE,add=TRUE)
    }
    )
  }
 
Tf=parms$Tf
for(t in 2:Tf){
	parms=list(lambda=c(run$t1[t,1],run$t2[t,1]),
		a=c(0,0),b=c(1,1)/4,alpha=c(1,1),
		beta=c(1,1)/5,M=150,X0=c())
	nullclines1(parms)
}

source("base-code.R")
parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=150,
  X0=c())
par(new=T)
nullclines(parms,main.text = "r=-1")
equi(parms,c(80,80))
equi(parms,c((parms$lambda[1]-1)/parms$alpha[1],0))
equi(parms,c(0,parms$lambda[2]-1)/parms$alpha[2])
equi(parms,c(10,80),col="white")
equi(parms,c(80,10),col="white")

points(run$x[,1],run$y[,1],type="b",pch=21,bg="gray",
       cex=0.5,col="gray")
#arrows(20,80,80,20,code=3,lwd=3,lty=1)
#legend("right",expression(italic(r)==-1),
#       cex=1.25, bty="n")
legend("topleft","A",bty="n",cex=1.5)



source("base-code.R")
cs<-seq(0,1,length=30)
parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=150,
  X0=c())

nullclines(parms,main.text = "r=0.95")

source("base-code-stochastic.R")
set.seed(3)
parms=list(b=c(1,1)/4,alpha=c(1,1),beta=c(1,1)/5,
           lambda=c(100,100),sd.log.lambda=c(0.05,0.05),
           rho=0.95,reps=2,Tf=55,a=c(0,0))
out=equi(parms,c(100,100))
run=um(parms=parms,x0=out$par)

Tf=parms$Tf
for(t in 2:Tf){
	parms=list(lambda=c(run$t1[t,1],run$t2[t,1]),
		a=c(0,0),b=c(1,1)/4,alpha=c(1,1),
		beta=c(1,1)/5,M=150,X0=c())
	nullclines1(parms)
}


source("base-code.R")
parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=150,
  X0=c())
par(new=T)
nullclines(parms)
equi(parms,c(80,80))
equi(parms,c((parms$lambda[1]-1)/parms$alpha[1],0))
equi(parms,c(0,parms$lambda[2]-1)/parms$alpha[2])
equi(parms,c(10,80),col="white")
equi(parms,c(80,10),col="white")

points(run$x[,1],run$y[,1],type="b",pch=21,bg="gray",
       cex=0.5,col="gray")
xs=out$par
#arrows(xs[1]*0.5,xs[2]*0.5,1.5*xs[1],1.5*xs[2],code=3,lwd=3,lty=1)
#legend("right",expression(italic(r)==0.95),
#       cex=1.25, bty="n")
legend("topleft","B",bty="n",cex=1.5)
dev.off()