source("base-code.R")
require(MASS)
require(matrixStats)
threshold=0.001 # extinction threshold
QE=TRUE # whether to use the traditional quasi-extinction threhold i.e. extinct first time crossing the threshold, or (FALSE) to only count as extinction if below threshold at the end. 
load.it=TRUE;load.file="FIGURE-stochasticbifs4-auto05.Rdata"

reps=200
mean.RI<-1.1/2

stable.equi=function(parms){
  with(as.list(parms),{
    f1=function(x,y)lambda[1]*x^2/(x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
    f2=function(x,y)lambda[2]*y^2/(y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
    L=10
    eq<-(lambda-1)/alpha
    c<-seq(0,1,length=L)
    xs<-c*eq[1]
    ys<-(1-c)*eq[2]
    Tf<-500
    for(t in 1:Tf){
      xs.new<-f1(xs,ys)
      ys.new<-f2(xs,ys)
      xs<-xs.new
      ys<-ys.new
    }
    if(length(which(xs*ys>0))>0){
      temp=min(which(xs*ys>0))
      return(cbind(xs[temp],ys[temp]))}
    else{return(NA)}
  }
  )
}


um=function(parms,x0){
  reps=parms$reps
  Tf=parms$Tf
  b=parms$b
  alpha=parms$alpha
  beta=parms$beta
  mean.lambda=parms$lambda
  sd.log.lambda=parms$sd.log.lambda
  rho=parms$rho
  rho.time=parms$rho.time
  log.mean.lambda=log(mean.lambda)-sd.log.lambda^2/2
  Sigma=diag(sd.log.lambda)
  Sigma[1,2]=rho*sd.log.lambda[1]
  Sigma[2,1]=rho*sd.log.lambda[1]
  x=matrix(x0[1],Tf,reps)
  y=matrix(x0[2],Tf,reps)
  t1=matrix(NA,Tf,reps)
  t2=matrix(NA,Tf,reps)
  z=x
  env=mvrnorm(n=reps,mu=c(0,0),Sigma=Sigma)
  for(t in 2:Tf){
   env=rho.time*env+sqrt(1-rho.time^2)*mvrnorm(n=reps,mu=c(0,0),Sigma=Sigma)
    temp=exp(cbind(env[,1]+log.mean.lambda[1],env[,2]+log.mean.lambda[2]))
    t1[t,]=temp[,1]
    t2[t,]=temp[,2]
    x[t,]=x[t-1,]^2*temp[,1]/(x[t-1,]+b[1]*y[t-1,])/(1+alpha[1]*x[t-1,]+beta[1]*y[t-1,])
    y[t,]=y[t-1,]^2*temp[,2]/(y[t-1,]+b[2]*x[t-1,])/(1+alpha[2]*y[t-1,]+beta[2]*x[t-1,])
    z[t,]=pmin(x[t,],y[t,])
  }
  return(list(x=x,y=y,z=z,t1=t1,t2=t2))
}


if(load.it){load(load.file)}else{
parms=list(
  lambda=c(100,100),
  b=c(0.1,0.1),
  alpha=c(1,1),
  beta=c(1,1)/10,
  M=120,reps=1000,Tf=100,sd.log.lambda=c(0.1,0.1),rho=0,rho.time=0) 

bs=seq(0.1,1,length=reps)
fitdiff<-seq(1,3,length=reps)
bifdatax<-matrix(0,reps,reps)
bifdatay<-matrix(0,reps,reps)
bifdataz<-matrix(0,reps,reps)

for(i in 1:reps){print(i)
  for(j in 1:reps){
    set.seed(1)
    tempoe=c(bs[i],0.15)
    parms$b=tempoe
    parms$lambda[1]=fitdiff[j]*parms$lambda[2]
    out=stable.equi(parms)
    tempx=NA
    tempy=NA
    tempz=NA
    if(length(out)>1){
      out2=um(parms,x0 = out)
      if(QE){tempx=length(which(colMins(out2$x)>threshold))/parms$reps
      tempy=length(which(colMins(out2$y)>threshold))/parms$reps
      tempz=length(which(colMins(out2$z)>threshold))/parms$reps}else{
        tempx=sum(sign(out2$x[parms$Tf,]))/parms$reps
        tempy=sum(sign(out2$y[parms$Tf,]))/parms$reps
        tempz=sum(sign(out2$z[parms$Tf,]))/parms$reps}
    }
    bifdatax[i,j]=tempx
    bifdatay[i,j]=tempy
    bifdataz[i,j]=tempz
  }
}
}

pdf("FIGURE-stochasticbifs4-auto05.pdf",width = 6,height=4)
par(mar=c(4,4.5,1,1),mfrow=c(1,1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(3,3,4,4,1,2,2,1),2,byrow=TRUE),heights=c(4,1.1))

##
## FIRST PAIR
##
par(mar=c(4,3,1,4),cex.lab=1.5,cex.axis=1.5)
plot(c(0,0),xlab="",ylab="",xaxt="n",yaxt="n",col="white",bty="n")
temp=seq(0,1,length=reps)
temp2=seq(1,0,length=reps)
image(temp,c(1,2),cbind(temp2,temp2),xaxt="n",yaxt="n",col=topo.colors(20),ylab="",xlab="",zlim=c(0,1))
axis(1,labels=TRUE)
mtext(text = "Probability of species loss",side=1,line=3,cex=1)




##
## SECOND TWO PANELS
##


par(mar=c(4,4.5,1,1),cex.lab=1.5,cex.axis=1.5)
image(bs/0.15,fitdiff,bifdatax,col=topo.colors(20),xlab=expression(paste("relative strength of PFD ",b[1]/b[2])),ylab=expression(paste("max fecundity ratio ",lambda[1]/lambda[2])),zlim=c(0,1))
legend("topleft","A Probability of extinction\n   of species 1",bty="n",cex=1.25)
image(bs/0.15,fitdiff,bifdatay,col=topo.colors(20),xlab=expression(paste("relative strength of PFD ",b[1]/b[2])),ylab=expression(paste("max fecundity ratio ",lambda[1]/lambda[2])),zlim=c(0,1))
legend("topleft","B Probability of extinction\n   of species 2",bty="n",cex=1.25)
dev.off()
#save.image("FIGURE-stochasticbifs4.Rdata")