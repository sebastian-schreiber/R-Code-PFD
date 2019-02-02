# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This code uses the deterministic and stochastic simulation code in base-code.R & base-code-stochastic.R
# To create Figure 3 in the main article and Figures S-1A,B in Supplement S3

# required code
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

pdf("why-rho-Supplement.pdf",width=8,height=4)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1), cex.lab = 1.5, cex.axis = 1.5)
nullclines(parms,main.text = expression(italic(r) == -1))
equi(parms,c(80,80))
equi(parms,c((parms$lambda[1]-1)/parms$alpha[1],0))
equi(parms,c(0,parms$lambda[2]-1)/parms$alpha[2])
equi(parms,c(10,80),col="white")
equi(parms,c(80,10),col="white")

source("base-code-stochastic.R")
set.seed(3)
parms=list(b=c(1,1)/4,alpha=c(1,1),beta=c(1,1)/5,
           lambda=c(100,100),sd.log.lambda=c(0.05,0.05),
           rho=-1,reps=2,Tf=100,a=c(0,0),rho.time=0.5,s=c(0,0))
out=equi(parms,c(100,100))
run=um(parms=parms,x0=out$par)
points(run$x[,1],run$y[,1],type="l",col="gray")
points(run$x[,1],run$y[,1],type="p",pch=21,bg="gray",
       cex=0.5,col="gray")
for(t in 1:parms$Tf){
	if(t %% 3 == 0){
		arrows(x0=run$x[t,1], y0=run$y[t,1], x1=run$x[t+1,1], y1=run$y[t+1,1], length = 0.1, angle = 20, col="gray")
	}
}
#arrows(20,80,80,20,code=3,lwd=3,lty=1)
#legend("right",expression(rho==-1),
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

nullclines(parms,main.text = expression(italic(r) == 0.95))
equi(parms,c(80,80))
equi(parms,c((parms$lambda[1]-1)/parms$alpha[1],0))
equi(parms,c(0,parms$lambda[2]-1)/parms$alpha[2])
equi(parms,c(10,80),col="white")
equi(parms,c(80,10),col="white")

source("base-code-stochastic.R")
set.seed(3)
parms=list(b=c(1,1)/4,alpha=c(1,1),beta=c(1,1)/5,
           lambda=c(100,100),sd.log.lambda=c(0.05,0.05),
           rho=0.95,reps=2,Tf=100,a=c(0,0),rho.time=0.5)
out=equi(parms,c(100,100))
run=um(parms=parms,x0=out$par)
#points(run$x[,1],run$y[,1],type="b",pch=21,bg="gray",
#       cex=0.5,col="gray")
points(run$x[,1],run$y[,1],type="l",col="gray")
points(run$x[,1],run$y[,1],type="p",pch=21,bg="gray",
       cex=0.5,col="gray")
for(t in 1:parms$Tf){
	if(t %% 3 == 0){
		arrows(x0=run$x[t,1], y0=run$y[t,1], x1=run$x[t+1,1], y1=run$y[t+1,1], length = 0.1, angle = 20, col="gray")
	}
}
xs=out$par
#arrows(xs[1]*0.5,xs[2]*0.5,1.5*xs[1],1.5*xs[2],code=3,lwd=3,lty=1)
#legend("right",expression(rho==0.95),
#       cex=1.25, bty="n")
legend("topleft","B",bty="n",cex=1.5)
dev.off()