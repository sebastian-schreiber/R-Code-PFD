# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This code uses the deterministic simulation code in base-code.R
# To create Figure 2 of the publication

# required code
source("base-code.R")
# create figure
reps<-200
asym.factor<-1/2
mean.RI<-1.1/2

cubic.condition=function(parms){
  with(parms,{
    ax=b[2]*beta[2]*lambda[1]
    bx=lambda[1]*(beta[2]+alpha[2]*b[2])-alpha[1]*lambda[2]
    cx=lambda[2]*(-beta[1]-alpha[1]*b[1])+alpha[2]*lambda[1]
    dx=-b[1]*beta[1]*lambda[2]
    disc=18*ax*bx*cx*dx-4*bx^3*dx+bx^2*cx^2-4*ax*cx^3-27*ax^2*dx^2
    return(disc)
    })

}


parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(0,0)+0.000001,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=120)
betas<-seq(0.00000001,1,length=reps)
fitdiff<-seq(1,3,length=reps)
bifdata<-matrix(0,reps,reps)
unstable.bifdata=bifdata
extra1<-matrix(0,reps,reps)
extra2<-matrix(0,reps,reps)
extra3<-matrix(0,reps,reps)
for(i in 1:reps){
  for(j in 1:reps){
    parms$lambda[1]=parms$lambda[2]*fitdiff[i]
    parms$beta=c(1,1)*betas[j]
    out=equi2(parms)
    temp=NA
    temp2=NA
    if(length(out)>1)temp=out[2]/(1+out[2])
    if(length(out)==1)temp2=out[1]/(1+out[1])
    unstable.bifdata[i,j]=temp2
    bifdata[i,j]=temp
 #   extra1[i,j]=(parms$lambda[1]/parms$lambda[2])-(parms$beta[1]+parms$b[1])/parms$alpha[1]
 #   extra2[i,j]=(parms$lambda[2]/parms$lambda[1])-(parms$beta[2]+parms$b[2])/parms$alpha[2]
    extra1[i,j]=(parms$lambda[1]/parms$alpha[1])*(parms$alpha[2]/parms$lambda[2])-
                parms$beta[1]/parms$alpha[1]-parms$b[1]
    extra2[i,j]=(parms$lambda[2]/parms$alpha[2])*(parms$alpha[1]/parms$lambda[1])-
                parms$beta[2]/parms$alpha[2]-parms$b[2] 
    extra3[i,j]=cubic.condition(parms)

  }
}

pdf("bifs3.pdf",width=7,height=7)
par(mar=c(4.5,4.5,1,1),mfrow=c(1,1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1,2,3,4,5,6),3,byrow=TRUE),heights=c(4,4,1.5))

image(betas,fitdiff,1-t(bifdata),col=topo.colors(20),ylab=expression(paste("max fecundity ratio ",
            lambda[1]/lambda[2])),xlab=expression(paste("strength of niche overlap")),zlim=c(0,1),main="No PFD")
#image(betas,fitdiff,1-t(unstable.bifdata),col=gray.colors(20),zlim=c(0,1),add=TRUE)
contour(betas,fitdiff,t(extra2),levels=c(0,0),add=TRUE,lty=2,lwd=2,drawlabels=FALSE)
contour(betas,fitdiff,t(extra3),levels=c(0,0),add=TRUE,lty=1,lwd=1,drawlabels=FALSE)
text(0.2,2,"coexistence",cex=1.5)
text(0.8,2,"exclusion",col="black",cex=1.5)
legend("topleft","A",bty="n",cex=1.5)
legend("topright",expression(paste("")),bty="n",cex=1.5)


parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/5,
  M=120)
betas<-seq(0.00000001,1,length=reps)
bs<-seq(0.00001,1,length=reps)
bifdata<-matrix(0,reps,reps)
extra1<-matrix(0,reps,reps)
extra2<-matrix(0,reps,reps)
for(i in 1:reps){
  for(j in 1:reps){
    parms$b=c(1,1)*bs[i]
    parms$beta=c(1,1)*betas[j]
    out=equi2(parms)
    temp=NA
    temp2=NA
    if(length(out)>1)temp=out[2]/(1+out[2])
    if(length(out)==1)temp2=out[1]/(1+out[1])
    bifdata[i,j]=temp
    unstable.bifdata[i,j]=temp2
    #extra1[i,j]=(parms$lambda[1]/parms$lambda[2])-(parms$beta[1]+parms$b[1])/parms$alpha[1]
    #extra2[i,j]=(parms$lambda[2]/parms$lambda[1])-(parms$beta[2]+parms$b[2])/parms$alpha[2]
    extra1[i,j]=(parms$lambda[1]/parms$alpha[1])*(parms$alpha[2]/parms$lambda[2])-
                parms$beta[1]/parms$alpha[1]-parms$b[1]
    extra2[i,j]=(parms$lambda[2]/parms$alpha[2])*(parms$alpha[1]/parms$lambda[1])-
                parms$beta[2]/parms$alpha[2]-parms$b[2]
    extra3[i,j]=cubic.condition(parms)

  }
}


image(bs,betas,1-t(bifdata),col=topo.colors(20),xlab=expression(paste("strength of PFD ",b[1]==b[2])),ylab=expression(paste("niche overlap ",rho)),zlim=c(0,1),main="Equal maximum fecundity ")
#image(bs,betas,1-t(unstable.bifdata),col=gray.colors(20),zlim=c(0,1),add=TRUE)
contour(bs,betas,extra1,levels=c(0,0),add=TRUE,lty=2,lwd=2,drawlabels=FALSE)
contour(bs,betas,extra2,levels=c(0,0),add=TRUE,lty=2,lwd=2,drawlabels=FALSE)
contour(bs,betas,extra3,levels=c(0,0),add=TRUE,lty=1,lwd=1,drawlabels=FALSE)
legend("topleft","B",bty="n",cex=1.5)
#legend("topright",expression(paste("no max fecundity difference")),bty="n",cex=1.5)
arrows(x0 = 0.34,y0 = 0.34,x1 = 0.5,y1 = 0.5,length=0.1,cex=1.5)
text(x = 0.4,y = 0.4,"non-additive\n\n contribution",srt=40,cex=0.9)

parms=list(
  lambda=c(100,100),
  a=c(0,0),
  b=c(0.1,1)/4,
  alpha=c(1,1),
  beta=c(1,1)/10,
  M=120)
bifdata<-matrix(0,reps,reps)

for(i in 1:reps){
  for(j in 1:reps){
    parms$b=c(1,1)*bs[i]
    parms$lambda[1]=fitdiff[j]*parms$lambda[2]
    out=equi2(parms)
    temp=NA
    temp2=NA
    if(length(out)>1)temp=out[2]/(1+out[2])
    if(length(out)==1)temp2=out[1]/(1+out[1])
    bifdata[i,j]=temp
    unstable.bifdata[i,j]=temp2
#    extra1[i,j]=(parms$lambda[1]/parms$lambda[2])-(parms$beta[1]+parms$b[1])/parms$alpha[1]
#    extra2[i,j]=(parms$lambda[2]/parms$lambda[1])-(parms$beta[2]+parms$b[2])/parms$alpha[2]
    extra1[i,j]=(parms$lambda[1]/parms$alpha[1])*(parms$alpha[2]/parms$lambda[2])-
                parms$beta[1]/parms$alpha[1]-parms$b[1]
    extra2[i,j]=(parms$lambda[2]/parms$alpha[2])*(parms$alpha[1]/parms$lambda[1])-
                parms$beta[2]/parms$alpha[2]-parms$b[2]
    extra3[i,j]=NA
    if(extra2[i,j]>0)extra3[i,j]=cubic.condition(parms)

  }
}

image(bs,fitdiff,1-bifdata,col=topo.colors(20),xlab=expression(paste("strength of PFD ",b[1]==b[2])),
  ylab=expression(paste("max fecundity ratio ",lambda[1]/lambda[2])),zlim=c(0,1),main="10% niche overlap, symmetric PFD")
#image(bs,fitdiff,1-(unstable.bifdata),col=gray.colors(20),zlim=c(0,1),add=TRUE)
#contour(bs,fitdiff,extra1,levels=c(0,0),add=TRUE,drawlabels=FALSE,lty=2,lwd=2)
contour(bs,fitdiff,extra2,levels=c(0,0),add=TRUE,drawlabels=FALSE,lty=2,lwd=2)
contour(bs,fitdiff,extra3,levels=c(0,0),add=TRUE,drawlabels=FALSE,lty=1,lwd=1,xlim=c(0,0.7))
#text(0.7,0.8*(max(fitdiff)-min(fitdiff))+min(fitdiff),"symmetric PFD",col="black",cex=1.5)
#text(0.7,0.7*(max(fitdiff)-min(fitdiff))+min(fitdiff),expression(lambda[2]==100),col="black",cex=1.5)
#text(0.7,0.9*(max(fitdiff)-min(fitdiff))+min(fitdiff),"10% niche overlap",col="black",cex=1.5)
legend("topleft","C",bty="n",cex=1.5)
#abline(v=mean.RI,lty=3,lwd=1)


bs=seq(0.1,1,length=reps)
for(i in 1:reps){
  for(j in 1:reps){
    tempoe=c(bs[i],0.15)
    parms$b=tempoe
    parms$lambda[1]=fitdiff[j]*parms$lambda[2]
    out=equi2(parms)
    temp=NA
    temp2=NA
    if(length(out)>1)temp=out[2]/(1+out[2])
    if(length(out)==1)temp2=out[1]/(1+out[1])
    bifdata[i,j]=temp
    unstable.bifdata[i,j]=temp2
     extra1[i,j]=(parms$lambda[1]/parms$alpha[1])*(parms$alpha[2]/parms$lambda[2])-
                parms$beta[1]/parms$alpha[1]-parms$b[1]
    extra2[i,j]=(parms$lambda[2]/parms$alpha[2])*(parms$alpha[1]/parms$lambda[1])-
                parms$beta[2]/parms$alpha[2]-parms$b[2]
    extra3[i,j]=cubic.condition(parms)
  }
}

image(bs/0.15,fitdiff,1-bifdata,col=topo.colors(20),xlab=expression(paste("relative PFD ",b[1]/b[2])),ylab=expression(paste("max fecunidty ratio ",lambda[1]/lambda[2])),zlim=c(0,1),main="10% niche overlap, asymmetric PFD")
#image(bs/0.15,fitdiff,1-(unstable.bifdata),col=gray.colors(20),zlim=c(0,1),add=TRUE)
contour(bs/0.15,fitdiff,extra1,levels=c(0,0),add=TRUE,drawlabels=FALSE,lty=2,lwd=2)
contour(bs/0.15,fitdiff,extra2,levels=c(0,0),add=TRUE,drawlabels=FALSE,lty=2,lwd=2)
contour(bs/0.15,fitdiff,extra3,levels=c(0,0),add=TRUE,drawlabels=FALSE,lty=1,lwd=1)
#text(5.2,0.8*(max(fitdiff)-min(fitdiff))+min(fitdiff),"asymmetric PFD",col="black",cex=1.5)
#text(5.2,0.7*(max(fitdiff)-min(fitdiff))+min(fitdiff),expression(b[2]==0.15),col="black",cex=1.5)
#text(5,0.9*(max(fitdiff)-min(fitdiff))+min(fitdiff),"10% niche overlap",col="black",cex=1.5)
legend("topleft","D",bty="n",cex=1.5)
#abline(v=1,lty=3,lwd=1)


par(mar=c(4,3,1,4))
temp=seq(0,1,length=reps)
temp2=seq(1,0,length=reps)
image(temp,c(1,2),cbind(temp2,temp2),xaxt="n",yaxt="n",col=topo.colors(20),ylab="",xlab="",zlim=c(0,1))
axis(1,labels=TRUE)
mtext(text = "freq. of species 1 at coexistence equilibrium",side=1,line=3,cex=1)

#image(temp,c(1,2),cbind(temp2,temp2),xaxt="n",yaxt="n",col=gray.colors(20),ylab="",xlab="",zlim=c(0,1))
#axis(1,labels=TRUE)
#mtext(text = "threshold frequency of species 1",side=1,line=3,cex=1)
dev.off()
