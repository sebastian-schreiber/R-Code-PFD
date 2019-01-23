# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This code uses the deterministic simulation code in base-code.R
# To create Figure 4 of the publication

# required code
source("base-code-stochastic-with-storage.R")
require(matrixStats)
# create figure
threshold=0.001 # extinction threshold
QE=TRUE # whether to use the traditional quasi-extinction threhold i.e. extinct first time crossing the threshold, or (FALSE) to only count as extinction if below threshold at the end. 

parms=list(b=c(1,1)/100,alpha=c(1,1),beta=c(1,1)/1,
           lambda=c(100,100),sd.log.lambda=c(0.5,0.5),rho=0,
           reps=1000,Tf=500,s=c(0.5,0.5),rho.time=0.5)
out=equi(parms,c(100,100))
x0=out$par


bif.slice=function(s){
  extinction=numeric(reps)
  igrs=numeric(reps)
  for(i in 1:reps){
    print(i)
    parms$rho=rhos[i]
    parms$rho.time=s
    #set.seed(1)
    run=um(parms=parms,x0=x0)
    if(QE){blerg=length(which(colMins(run$z)<0.001))/parms$reps}else{blerg=length(which(run$z[parms$Tf,]<0.001))/parms$reps} 
    extinction[i]=blerg
    igrs[i]=igr(parms)
  }
  return(list(extinction=extinction,igr=igrs))
}

reps=17 #17
rhos=seq(-0.9,0.9,length=reps)
parms$reps=20000 # use 20,000 for final figure
parms$Tf=50 # 50 in main text figure 
parms$sd.log.lambda=c(1.5,1.5)
reps3=10 #10
storage.stuff=matrix(NA,reps,reps3)
igr.stuff=matrix(NA,reps,reps3)
ss=seq(0,0.9,length=reps3)
for(i in 1:reps3){
  print(i)
  print(i)
  out=bif.slice(ss[i])
  storage.stuff[,i]=1-out$extinction
  igr.stuff[,i]=out$igr}


#save.image("FIGURE-stochastic-bif-with-storage.Rdata")
load("FIGURE-stochastic-bif-with-storage.Rdata")
pdf("FIGURE-stochastic-with-storage.pdf",width = 7,height=4)
layout(mat = cbind(c(1),c(2),c(3)),widths = c(3.5,3.5,1.25))
par(mar=c(4,4,1,1),cex.axis=1.25,cex.lab=1.25)

matplot(rhos,igr.stuff,type="l",lty=1,lwd=3,
        bty="n",xlab="cross-correlation",
        ylab="invasion growth rate",
        col=topo.colors(reps3))
legend("topleft","A",bty="n",cex=1.25)
arrows(x0=-0,y0=min(igr.stuff),x1=-0,y1=max(igr.stuff),length = 0.2,lwd=2,lty=1,col=gray(0.3))
text(-0.,0.2,"increasing\n autocorrelation",srt=90,cex=1,col=gray(0.3))
matplot(rhos,storage.stuff,type="l",ylim=c(0,1),lty=1,lwd=3,
        bty="n",xlab="cross-correlation",
        ylab="persistence probability",
        col=topo.colors(reps3))
legend("topleft","B",bty="n",cex=1.25)
arrows(x0=-0,y0=1,x1=-0,y1=0,length = 0.2,lwd=2,lty=1,col=gray(0.3))
text(-0.,0.5,"increasing\n autocorrelation",srt=90,cex=1,col=gray(0.3))
par(mar=c(4,1,1,4))
temp=ss
temp2=ss
image(c(1,2),temp,rbind(temp2,temp2),xaxt="n",yaxt="n",col=topo.colors(reps3),ylab="",xlab="",zlim=c(0,1))
axis(4,labels=TRUE)
mtext(text = "temporal autocorrelation",side=4,line=3,cex=0.75)
dev.off()