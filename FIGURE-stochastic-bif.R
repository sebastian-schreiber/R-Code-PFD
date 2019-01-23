# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This code uses the deterministic simulation code in base-code.R
# To create Figure 5 of the publication

# required code  
source("base-code-stochastic.R")
library(matrixStats)
# create figure  
  threshold=0.001 # extinction threshold
  QE=TRUE # whether to use the traditional quasi-extinction threhold i.e. extinct first time crossing the threshold, or (FALSE) to only count as extinction if below threshold at the end. 
  load.it=TRUE;load.file="FIGURE-stochastic-bifauto05.Rdata" #files FIGURE-stochastic-bifauto05.Rdata and FIGURE-stochastic-bif.Rdata
  
  parms=list(b=c(1,1)/4,alpha=c(1,1),beta=c(1,1)/5,
             lambda=c(100,100),sd.log.lambda=c(0.5,0.5),rho=0,reps=5,Tf=100,a=c(0,0),rho.time=0)
  out=equi(parms,c(100,100))
  x0=out$par
  reps=40 #40 in final
  rhos<-seq(-0.95,0.95,length=reps)
  sigmas<-seq(0,1,length=reps)
  
  
  bif.function=function(b=c(1,1)/4){
    parms$b=b
    extinction=matrix(NA,reps,reps)
    guess=c(100,100)
    for(i in 1:reps){
      for(j in 1:reps){
        parms$Tf=50
        parms$sd.log.lambda=sigmas[i]*c(1,1)
        parms$rho=rhos[j]
        out=stable.equi(parms)
        x0=out
        #set.seed(1)
        run=um(parms=parms,x0=x0)
        if(QE){blerg=length(which(colMins(run$z)<threshold))/parms$reps}else{
        blerg=length(which(run$z[parms$Tf,]<threshold))/parms$reps}
        extinction[i,j]=blerg
      }
    }
    return(extinction)
  }
  

  pdf("FIGURE-stochastic-bif.pdf",width = 7,height=6)
  par(mar=c(4.5,4.5,1,1),cex.lab=1.5,cex.axis=1.25)
  layout(matrix(c(1,3,2,3,0,4),3,byrow=TRUE),heights=c(4,4,1.5),widths=c(2,4))
  
  out=equi(parms,c(100,100))
  x0=out$par
  parms$rho=0.6
  parms$Tf=100
  parms$reps=10000 # final version 10000
  set.seed(1);out=um(parms,x0)
  matplot(out$z,type="l",lty=1,xlab="years",ylab="min. density")
  legend("topright",expression(r==0.6),bty="n")
  legend("topleft","A",bty="n",cex=1.5)
  parms$rho=-0.6
  set.seed(1);out=um(parms,x0)
  matplot(out$z,type="l",lty=1,xlab="years",ylab="min. density")
  legend("topright",expression(r==-0.6),bty="n")
  legend("topleft","B",bty="n",cex=1.5)
  
  if(load.it){load(load.file)}else{
  extinction=bif.function()}
  zz=1-extinction
  #save.image("FIGURE-stochastic-bifauto05.Rdata")
  image(sigmas,rhos,zz,col=topo.colors(20),xlab=expression(paste("environmental variance ",sigma^2)),ylab=expression(paste("cross-correlation ",r)),zlim=c(0,1))
  legend("topleft","C",bty="n",cex=1.5)
  
  
  par(mar=c(4,3,1,4))
  temp=seq(0,1,length=reps)
  temp2=seq(1,0,length=reps)
  image(temp,c(1,2),cbind(temp2,temp2),xaxt="n",yaxt="n",col=topo.colors(20),ylab="",xlab="",zlim=c(0,1))
  axis(1,labels=TRUE)
  mtext(text = "probability of species loss",side=1,line=3,cex=1)
  dev.off()