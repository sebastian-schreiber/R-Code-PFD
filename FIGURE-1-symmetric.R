# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This code uses the deterministic simulation code in base-code.R
# To create Figure 1 of the publication

# required code
source("base-code.R")

# create figure
cs<-seq(0,1,length=30)
parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=120,
  X0=cbind(cs*1,(1-cs)*1),
  Xarrows=seq(2,30,by=3))

pdf("symmetric.pdf",width=7,height=7)

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 1, 1), cex.lab = 1.5, cex.axis = 1.5)
logfit(parms,approx = FALSE) # To see approximations switch FALSE to TRUE
legend("topleft","A",bty="n",cex=1.5)
nullclines(parms)
legend("topleft","B",bty="n",cex=1.5)
equi(parms,c(80,80))
equi(parms,c((parms$lambda[1]-1)/parms$alpha[1],0))
equi(parms,c(0,parms$lambda[2]-1)/parms$alpha[2])
equi(parms,c(10,80),col="white")
equi(parms,c(80,10),col="white")

parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(3,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=120,
  X0=cbind(cs*1,(1-cs)*1),
  Xarrows=seq(2,30,by=3))
logfit(parms,leg=FALSE,approx = FALSE)
legend("topleft","C",bty="n",cex=1.5)
nullclines(parms)
legend("topleft","D",bty="n",cex=1.5)
equi(parms,c(80,30),col="white")
equi(parms,c((parms$lambda[1]-1)/parms$alpha[1],0),col="black")
equi(parms,c(0,parms$lambda[2]-1)/parms$alpha[2],col="black")
dev.off()
