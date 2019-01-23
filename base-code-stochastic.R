# This R code is for the paper 
# "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity"
# by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss
# that has been accepted for publication in Ecology. 

# This file providews the  base commands 
# for simulatting and exploring the stochastic, two species competition model model  
# Ni[t+1]=Ni[t]*(Ni[t]/(a[i]+Ni[t]+bi*Nj[t]))*(lambdai[t]/(1+alphai*Ni[t]+betai*Nj[t]))
# where lambdai is the intrinsic rate of growth of species i
# alphai/betai are the strength of intra/interspecific competition on species i
# ai and bi determine the strength of positive frequency dependence
# Most of this code (as with the Ecology paper) only exame the model with ai=0


# required packages
require(rootSolve)
require(MASS)


# The following code solves for 
# equilibria given an initial guess
# input: list of parms lambda, a, b, alpha, beta and initial guess
# output: an equilibrium

equi=function(parms,guess){
  with(as.list(parms),{
    f1=function(x,y)lambda[1]*x^2/(a[1]+x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
    f2=function(x,y)lambda[2]*y^2/(a[2]+y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
    h=function(x){(f1(x[1],x[2])-x[1])^2+(f2(x[1],x[2])-x[2])^2}
    out=optim(guess,h)
    return(out)
  }
  )
}

# The following numerically approximates stable equilibria when they exist
# input: list of parms lambda, a, b, alpha, beta
# output: an equilibrium when it exists else NA
stable.equi=function(parms){
  with(as.list(parms),{
    #assuming that a[i]=0
    f1=function(x,y)lambda[1]*x^2/(a[1]+x+b[1]*y)/(1+alpha[1]*x+beta[1]*y)
    f2=function(x,y)lambda[2]*y^2/(a[2]+y+b[2]*x)/(1+alpha[2]*y+beta[2]*x)
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


# this function simulates the stochastic model
# inputs: parameters include length Tf of simulation, mean value of lambdas, standard deviation of 
# log lambdas, autocorrelation rho.time, etc. and initial condition x0
# output: simulation as list of species densities, minimum z of species densities, and lambda values (t1 and t2) 
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

