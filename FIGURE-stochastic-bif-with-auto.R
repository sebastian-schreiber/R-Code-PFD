source("base-code-stochastic-with-storage.R")
require(matrixStats)
threshold=0.001 # extinction threshold
QE=TRUE # whether to use the traditional quasi-extinction threhold i.e. extinct first time crossing the threshold, or (FALSE) to only count as extinction if below threshold at the end. 

parms=list(b=c(1,1)/10,alpha=c(1,1),beta=c(1,1)/5,
           lambda=c(100,100),sd.log.lambda=c(0.2,0.2),rho=0,
           reps=1000,Tf=500,s=c(0,0),rho.time=0.5)
out=equi(parms,c(100,100))
x0=out$par


bif.slice=function(s){
  extinction=numeric(reps)
  igrs=numeric(reps)
  for(i in 1:reps){
    print(i)
    parms$rho=rhos[i]
    parms$rho.time=s
    run=um(parms=parms,x0=x0)
    if(QE){blerg=length(which(colMins(run$z)<0.001))/parms$reps}else{blerg=length(which(run$z[parms$Tf,]<0.001))/parms$reps}
    extinction[i]=blerg
  }
  return(list(extinction=extinction))
}

reps=17 #17
rhos=seq(-0.9,0.9,length=reps)
parms$reps=20000 # use 20,000 for final figure
parms$Tf=50 # 50 in main text figure, 25 for low storage figure 
parms$sd.log.lambda=c(0.5,0.5)
reps3=10 #10
storage.stuff=matrix(NA,reps,reps3)
igr.stuff=matrix(NA,reps,reps3)
rho.times=seq(0,0.9,length=reps3)
for(i in 1:reps3){
  print(i)
  out=bif.slice(rho.times[i])
  storage.stuff[,i]=1-out$extinction
  }


#save.image("FIGURE-stochastic-bif-with-auto.Rdata")
#load("FIGURE-stochastic-bif-with-auto.Rdata")
pdf("FIGURE-stochastic-bif-with-auto.pdf",width = 6,height=4)
layout(mat = cbind(c(1),c(2)),widths = c(3.5,0.8))
par(mar=c(4,4,1,1),cex.axis=1.25,cex.lab=1.25)
matplot(rhos,storage.stuff,type="l",ylim=c(0,max(storage.stuff)),lty=1,lwd=3,
        bty="n",xlab="cross-correlation",
        ylab="persistence probability",
        col=topo.colors(reps3))
par(mar=c(4,1,1,4))
legend("topleft","C",bty="n",cex=1.25)
temp=rho.times
temp2=rho.times
image(c(1,2),temp,rbind(temp2,temp2),xaxt="n",yaxt="n",col=topo.colors(reps3),ylab="",xlab="",zlim=c(0,1))
axis(4,labels=TRUE)
mtext(text = "temporal autocorrelation",side=4,line=3,cex=1)
dev.off()